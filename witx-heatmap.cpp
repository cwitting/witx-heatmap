// #include <omp.h>

#include <valhalla/config.h>
#include <valhalla/tyr/actor.h>

#include <boost/property_tree/ptree.hpp>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

#include "json.hpp"

struct Coordinate {
  double lat{};
  double lon{};
};

using Route = std::vector<Coordinate>;

struct MatchedRoute {
  Route route;
  std::vector<std::string> osm_route;  // Store the OSM route as a way ids
};

#define ROUTE_FILE "/home/christian/git/witx-heatmap/data/christian_strava.geojson"
#define VALHALLA_CONFIG_FILE "/home/christian/git/witx-heatmap/data/routing/valhalla_data/valhalla.json"

std::vector<Route> load_all_routes() {
  std::vector<Route> routes;

  std::ifstream file(ROUTE_FILE);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open route file: " ROUTE_FILE);
  }

  nlohmann::json doc;
  file >> doc;

  for (const auto& feature : doc["features"]) {
    const auto& geometry = feature["geometry"];
    if (geometry["type"].get<std::string>() != "LineString") {
      continue;
    }

    Route route;
    for (const auto& point : geometry["coordinates"]) {
      // GeoJSON coordinates are [lon, lat]
      route.push_back(Coordinate{point[1].get<double>(), point[0].get<double>()});
    }
    routes.push_back(std::move(route));
  }

  return routes;
}
// Using chrono timer
class TimerLog {
 public:
  TimerLog(const std::string& name) : name_(name), start_(std::chrono::high_resolution_clock::now()) {}

  ~TimerLog() { stop(); }

  void stop() {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start_;
    fprintf(stderr, "%s took %.6f seconds\n", name_.c_str(), elapsed.count());
  }

 private:
  std::string name_;
  std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

class RouteMatcher {
 public:
  RouteMatcher() {
    fprintf(stderr, "Initializing Valhalla with config from: %s\n", VALHALLA_CONFIG_FILE);
    TimerLog timer("Valhalla initialization");
    const auto& config = valhalla::config(VALHALLA_CONFIG_FILE);
    // auto_cleanup releases the loki/thor/odin workers' caches between calls.
    actor_ = std::make_unique<valhalla::tyr::actor_t>(config, /*auto_cleanup=*/true);
  }

  // Match routes to the road network using Valhalla's trace_route (map matching)
  std::optional<MatchedRoute> matchRoute(const Route& route) {
    nlohmann::json request;
    request["costing"] = "pedestrian";
    for (const auto& coord : route) {
      request["shape"].push_back({{"lat", coord.lat}, {"lon", coord.lon}});
    }

    valhalla::Api api;
    std::string json_str;
    try {
      json_str = actor_->trace_route(request.dump(), nullptr, &api);
    } catch (const std::exception& e) {
      fprintf(stderr, "Valhalla trace_route failed: %s\n", e.what());
      return std::nullopt;
    }

    // print the raw json result to stderr
    fprintf(stderr, "Valhalla trace_route result: %s\n", json_str.c_str());

    MatchedRoute matched_route;
    matched_route.route = route;

    // Add the OSM way ids to the matched_route.osm_route vector
    for (const auto& trip_route : api.trip().routes()) {
      for (const auto& leg : trip_route.legs()) {
        for (const auto& node : leg.node()) {
          if (node.edge().way_id() != 0) {
            matched_route.osm_route.push_back(std::to_string(node.edge().way_id()));
          }
        }
      }
    }

    return matched_route;
  }

  std::vector<MatchedRoute> matchAllRoutes(const std::vector<Route>& routes) {
    std::vector<MatchedRoute> matched_routes;
    int matched_count = 0;
    int total_count = 0;
    for (const auto& route : routes) {
      TimerLog timer("Matching route " + std::to_string(matched_count + 1) + "/" + std::to_string(total_count + 1) +
                     "/" + std::to_string(routes.size()));
      auto matched_route = matchRoute(route);
      total_count++;
      if (!matched_route) {
        continue;
      }
      matched_routes.push_back(*matched_route);
      matched_count++;
      if (matched_count >= 30) {
        break;
      }
    }
    return matched_routes;
  }

 private:
  std::unique_ptr<valhalla::tyr::actor_t> actor_;

  // Decode polyline geometry
  std::vector<Coordinate> decodePolyline(const std::string& encoded) const {
    std::vector<Coordinate> coordinates;
    int index = 0;
    int len = encoded.length();
    int lat = 0;
    int lng = 0;

    while (index < len) {
      int b;
      int shift = 0;
      int result = 0;

      do {
        b = encoded[index++] - 63;
        result |= (b & 0x1f) << shift;
        shift += 5;
      } while (b >= 0x20);

      int dlat = ((result & 1) ? ~(result >> 1) : (result >> 1));
      lat += dlat;

      shift = 0;
      result = 0;

      do {
        b = encoded[index++] - 63;
        result |= (b & 0x1f) << shift;
        shift += 5;
      } while (b >= 0x20);

      int dlng = ((result & 1) ? ~(result >> 1) : (result >> 1));
      lng += dlng;

      coordinates.push_back(Coordinate{lat / 1e5, lng / 1e5});
    }

    return coordinates;
  }
};

// ============================================================================
// Main program
// ============================================================================
int main(int argc, char* argv[]) {
  // Print program name
  fprintf(stderr, "Witx Heatmap Route Planner\n");
  auto routes = load_all_routes();
  RouteMatcher matcher;
  auto matched_routes = matcher.matchAllRoutes(routes);
  fprintf(stderr, "Matched %zu routes out of %zu\n", matched_routes.size(), routes.size());
  for (const auto& matched_route : matched_routes) {
    fprintf(stderr, "Matched route with %zu way ids\n", matched_route.osm_route.size());
    for (const auto& way_id : matched_route.osm_route) {
      std::cout << way_id << " ";
    }
  }
  std::cout << std::endl;
}
