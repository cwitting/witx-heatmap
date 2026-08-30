// #include <omp.h>

#include <valhalla/config.h>
#include <valhalla/tyr/actor.h>

#include <boost/property_tree/ptree.hpp>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <opencv4/opencv2/opencv.hpp>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "httplib.h"
#include "json.hpp"

struct Coordinate {
  double lat{};
  double lon{};
};

using Route = std::vector<Coordinate>;

struct WaySegment {
  std::string way_id;
  uint64_t edge_id{};  // Valhalla's internal directed-edge id: stable identity for the exact
                       // maximal way-segment traversed, safe to use as a traversal-count key
  Route geometry;      // The portion of the way's geometry actually traversed by the route
};

struct MatchedRoute {
  Route route;
  std::vector<WaySegment> way_segments;
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
    request["shape_match"] = "map_snap";
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
    // fprintf(stderr, "Valhalla trace_route result: %s\n", json_str.c_str());

    MatchedRoute matched_route;
    matched_route.route = route;

    // For each traversed edge, slice out the sub-range of the leg shape it covers so we
    // know exactly which part of the way (not just which way) was used.
    for (const auto& trip_route : api.trip().routes()) {
      for (const auto& leg : trip_route.legs()) {
        auto leg_shape = decodePolyline(leg.shape());
        for (const auto& node : leg.node()) {
          const auto& edge = node.edge();
          if (edge.way_id() == 0 || leg_shape.empty()) {
            continue;
          }
          size_t begin = std::min<size_t>(edge.begin_shape_index(), leg_shape.size() - 1);
          size_t end = std::min<size_t>(edge.end_shape_index(), leg_shape.size() - 1);
          Route geometry(leg_shape.begin() + begin, leg_shape.begin() + end + 1);
          matched_route.way_segments.push_back({std::to_string(edge.way_id()), edge.id(), std::move(geometry)});
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

      // Valhalla's default shape_format is polyline6
      coordinates.push_back(Coordinate{lat / 1e6, lng / 1e6});
    }

    return coordinates;
  }
};

cv::Mat generate_placeholder_tile(int z, int x, int y) {
  // Generate a simple 256x256 PNG image with a solid color
  const int width = 256;
  const int height = 256;
  cv::Mat png_data(height, width, CV_8UC4, cv::Scalar(255, 255, 255, 255));

  // Fill with a color based on the tile coordinates for demonstration
  unsigned char r = static_cast<unsigned char>((x * 37) % 256);
  unsigned char g = static_cast<unsigned char>((y * 59) % 256);
  unsigned char b = static_cast<unsigned char>((z * 83) % 256);

  for (int y_idx = 0; y_idx < height; ++y_idx) {
    for (int x_idx = 0; x_idx < width; ++x_idx) {
      cv::Vec4b& pixel = png_data.at<cv::Vec4b>(y_idx, x_idx);
      pixel[0] = b;    // Blue
      pixel[1] = g;    // Green
      pixel[2] = r;    // Red
      pixel[3] = 255;  // Alpha
    }
  }

  return png_data;
}

int main(int argc, char** argv) {
  // Create HTTP server
  httplib::Server svr;

  // Add CORS middleware to all responses
  svr.set_post_routing_handler([](const httplib::Request&, httplib::Response& res) {
    res.set_header("Access-Control-Allow-Origin", "*");
    res.set_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
    res.set_header("Access-Control-Allow-Headers", "Content-Type");
  });

  // Handle OPTIONS requests for CORS preflight
  svr.Options(".*", [](const httplib::Request&, httplib::Response& res) {
    res.set_header("Access-Control-Allow-Origin", "*");
    res.set_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
    res.set_header("Access-Control-Allow-Headers", "Content-Type");
    res.status = 204;
  });

  // Health check endpoint
  svr.Get("/health", [](const httplib::Request&, httplib::Response& res) {
    res.set_content("{\"status\": \"ok\"}", "application/json");
  });

  std::cout << "Loading OSM data from data/zealand.pbf..." << std::endl;

  // Print program name
  fprintf(stderr, "Witx Heatmap Route Planner\n");
  auto routes = load_all_routes();
  RouteMatcher matcher;
  auto matched_routes = matcher.matchAllRoutes(routes);

  // Main route planning endpoint
  svr.Get(R"(/tiles/(\d+)/(\d+)/(\d+).png)", [&matched_routes](const httplib::Request& req, httplib::Response& res) {
    // Heatmap XYZ tile request from url like /tiles/{z}/{x}/{y}.png
    for (const auto& param : req.path_params) {
      fprintf(stderr, "Path param: %s = %s\n", param.first.c_str(), param.second.c_str());
    }
    int z = std::stoi(req.matches[1]);
    int x = std::stoi(req.matches[2]);
    int y = std::stoi(req.matches[3]);
    fprintf(stderr, "Received tile request for z=%d, x=%d, y=%d\n", z, x, y);

    // For now just return a placeholder PNG image (256x256 pixel)
    cv::Mat png_data = generate_placeholder_tile(z, x, y);
    std::vector<unsigned char> buffer;
    cv::imencode(".png", png_data, buffer);
    res.set_content(reinterpret_cast<const char*>(buffer.data()), buffer.size(), "image/png");
  });

  svr.set_exception_handler([](const auto& req, auto& res, std::exception_ptr ep) {
    auto fmt = "<h1>Error 500</h1><p>%s</p>";
    char buf[BUFSIZ];
    fprintf(stderr, "Exception caught handling request %s\n", req.path.c_str());
    try {
      std::rethrow_exception(ep);
      fprintf(stderr, "A Exception handling request %s\n", req.path.c_str());
    } catch (std::exception& e) {
      snprintf(buf, sizeof(buf), fmt, e.what());
      fprintf(stderr, "Exception handling request %s: %s\n", req.path.c_str(), e.what());
    } catch (...) {  // See the following NOTE
      snprintf(buf, sizeof(buf), fmt, "Unknown Exception");
      fprintf(stderr, "Unknown exception handling request %s\n", req.path.c_str());
    }
    res.set_content(buf, "text/html");
    res.status = httplib::StatusCode::InternalServerError_500;
  });

  // Start server
  std::cout << "\n==================================================" << std::endl;
  std::cout << "Route Planning API Server" << std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << "Server starting on http://0.0f.0f.0f:8080" << std::endl;
  std::cout << "\nEndpoints:" << std::endl;
  std::cout << "  GET /health - Health check" << std::endl;
  std::cout << "  GET /plan?lat=55.59813&lon=11.97320&depth=10&radius=10&samples=8&target_distance=50" << std::endl;
  std::cout << "\nParameters:" << std::endl;
  std::cout << "  lat              - Origin latitude (default: 55.59813)" << std::endl;
  std::cout << "  lon              - Origin longitude (default: 11.97320)" << std::endl;
  std::cout << "  depth            - Tree depth (default: 10)" << std::endl;
  std::cout << "  radius           - Search radius in km (default: 10)" << std::endl;
  std::cout << "  samples          - Samples per node (default: 8)" << std::endl;
  std::cout << "  target_distance  - Target route distance in km (default: 50)" << std::endl;
  std::cout << "  distance_deviation - Distance deviation tolerance (default: 5000)" << std::endl;
  std::cout << "  max_overlap      - Max overlap percentage (default: 0.0f2)" << std::endl;
  std::cout << "  min_gravel       - Min gravel percentage (default: 0.2)" << std::endl;
  std::cout << "  format           - Output format: geojson or gpx (default: geojson)" << std::endl;
  std::cout << "==================================================" << std::endl;

  svr.listen("0.0.0.0", 8080);

  return 0;
}

// ============================================================================
// Main program
// ============================================================================
// int main(int argc, char* argv[]) {
//   // Print program name
//   fprintf(stderr, "Witx Heatmap Route Planner\n");
//   auto routes = load_all_routes();
//   RouteMatcher matcher;
//   auto matched_routes = matcher.matchAllRoutes(routes);
//   // OsmFinder osm_finder("/home/christian/git/witx-heatmap/data/routing/valhalla_data/merged.osm.pbf");

//   fprintf(stderr, "Matched %zu routes out of %zu\n", matched_routes.size(), routes.size());

//   // Count how many times each physical way-segment (Valhalla directed edge) was traversed
//   // across all matched routes. edge_id is stable across matches as long as the tiles don't
//   // change, and distinguishes direction, unlike way_id which can be shared by many segments.
//   std::unordered_map<uint64_t, int> traversal_counts;
//   for (const auto& matched_route : matched_routes) {
//     for (const auto& segment : matched_route.way_segments) {
//       traversal_counts[segment.edge_id]++;
//     }
//   }

//   for (const auto& matched_route : matched_routes) {
//     fprintf(stderr, "Matched route with %zu way segments\n", matched_route.way_segments.size());
//     for (const auto& segment : matched_route.way_segments) {
//       fprintf(stderr, "Way %s (edge %llu, traversed %d times): %zu coordinates\n", segment.way_id.c_str(),
//               static_cast<unsigned long long>(segment.edge_id), traversal_counts[segment.edge_id],
//               segment.geometry.size());
//       for (const auto& coord : segment.geometry) {
//         std::cerr << coord.lat << "," << coord.lon << " ";
//       }
//       std::cerr << std::endl;
//     }
//     std::cerr << std::endl;
//   }
//   std::cerr << std::endl;
// }
