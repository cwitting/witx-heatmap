// #include <omp.h>

#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <osrm/coordinate.hpp>
#include <osrm/engine_config.hpp>
#include <osrm/json_container.hpp>
#include <osrm/nearest_parameters.hpp>
#include <osrm/osrm.hpp>
#include <osrm/route_parameters.hpp>
#include <osrm/status.hpp>
#include <stdexcept>
#include <vector>

#include "engine/api/match_parameters.hpp"
#include "json.hpp"
#include "util/json_renderer.hpp"

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
#define OSRM_DATA_DIR "/home/christian/git/witx-heatmap/data/routing/data/merged"

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
    // Configure OSRM
    osrm::EngineConfig config;

    config.storage_config = {OSRM_DATA_DIR};
    config.use_shared_memory = false;
    config.algorithm = osrm::EngineConfig::Algorithm::MLD;

    // Create OSRM instance
    fprintf(stderr, "Initializing OSRM with data from: %s\n", OSRM_DATA_DIR);
    TimerLog timer("OSRM initialization");
    osrm_ = std::make_unique<osrm::OSRM>(config);
  }

  // Match routes to the road network using OSRM
  std::optional<MatchedRoute> matchRoute(const Route& route) {
    osrm::RouteParameters params;
    for (const auto& coord : route) {
      params.coordinates.emplace_back(osrm::util::FloatLongitude(coord.lon), osrm::util::FloatLatitude(coord.lat));
    }

    osrm::json::Object result;
    osrm::MatchParameters match_params;
    match_params.coordinates = params.coordinates;
    match_params.annotations = true;
    match_params.annotations_type = osrm::MatchParameters::AnnotationsType::Nodes;
    const auto status = osrm_->Match(match_params, result);

    if (status != osrm::Status::Ok) {
      fprintf(stderr, "OSRM Match failed with status: %d\n", static_cast<int>(status));
      return std::nullopt;
    }

    MatchedRoute matched_route;
    matched_route.route = route;

    // print the raw json result to stderr
    std::string json_str;
    osrm::util::json::render(json_str, result);
    fprintf(stderr, "OSRM Match result: %s\n", json_str.c_str());

    // Add the osm way ids to the matched_route.osm_route vector

    if (result.values.count("matchings") > 0) {
      const auto& matchings = std::get<osrm::util::json::Array>(result.values["matchings"]);
      for (const auto& matching : matchings.values) {
        const auto& matching_obj = std::get<osrm::util::json::Object>(matching);
        if (matching_obj.values.count("legs") > 0) {
          const auto& legs = std::get<osrm::util::json::Array>(matching_obj.values.at("legs"));
          for (const auto& leg : legs.values) {
            const auto& leg_obj = std::get<osrm::util::json::Object>(leg);
            if (leg_obj.values.count("annotation") > 0) {
              const auto& annotation = std::get<osrm::util::json::Object>(leg_obj.values.at("annotation"));
              if (annotation.values.count("nodes") > 0) {
                const auto& nodes = std::get<osrm::util::json::Array>(annotation.values.at("nodes"));
                for (const auto& node : nodes.values) {
                  matched_route.osm_route.push_back(std::to_string(std::get<osrm::util::json::Number>(node).value));
                }
              }
            }
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
  std::unique_ptr<osrm::OSRM> osrm_;

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
  // for (const auto& matched_route : matched_routes) {
  //   fprintf(stderr, "Matched route with %zu way ids\n", matched_route.osm_route.size());
  //   for (const auto& way_id : matched_route.osm_route) {
  //     std::cout << way_id << " ";
  //   }
  // }
  std::cout << std::endl;
}
