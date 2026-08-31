// #include <omp.h>

#include <valhalla/config.h>
#include <valhalla/tyr/actor.h>

#include <boost/property_tree/ptree.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <opencv2/core/mat.hpp>
#include <opencv4/opencv2/opencv.hpp>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "httplib.h"
#include "json.hpp"

#define ROUTE_FILE "/home/christian/git/witx-heatmap/data/christian_strava2.geojson"
#define VALHALLA_CONFIG_FILE "/home/christian/git/witx-heatmap/data/routing/valhalla_data/valhalla.json"
#define HEATMAP_FILE "/home/christian/git/witx-heatmap/data/heatmap.json"
#define RESTORE 0

struct Coordinate {
  double lat{};
  double lon{};
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Coordinate, lat, lon)

using Route = std::vector<Coordinate>;

struct WaySegment {
  std::string way_id;
  uint64_t edge_id{};     // Valhalla's internal directed-edge id: stable identity for the exact
                          // maximal way-segment traversed, safe to use as a traversal-count key
  Route geometry;         // The portion of the way's geometry actually traversed by the route
  int traversal_count{};  // How many times this way-segment was traversed across all matched routes

  bool inBBox(double min_lat, double min_lon, double max_lat, double max_lon) const {
    for (int i = 0; i < geometry.size(); ++i) {
      const auto& coord = geometry[i];
      if (coord.lat >= min_lat && coord.lat <= max_lat && coord.lon >= min_lon && coord.lon <= max_lon) {
        return true;
      }
      if (i > 0) {
        const auto& prev_coord = geometry[i - 1];
        // Check if the line segment between prev_coord and coord intersects the bounding box
        if ((prev_coord.lat < min_lat && coord.lat > max_lat) || (prev_coord.lat > max_lat && coord.lat < min_lat) ||
            (prev_coord.lon < min_lon && coord.lon > max_lon) || (prev_coord.lon > max_lon && coord.lon < min_lon)) {
          return true;
        }
      }
    }
    return false;
  }
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(WaySegment, way_id, edge_id, geometry, traversal_count)

std::pair<double, double> latLonToMeters(double lat, double lon) {
  double x = lon * 20037508.34 / 180.0;
  double y = std::log(std::tan((90.0 + lat) * M_PI / 360.0)) / (M_PI / 180.0);
  y = y * 20037508.34 / 180.0;
  return {x, y};
}

std::pair<double, double> metersToLatLon(double x, double y) {
  double lon = (x / 20037508.34) * 180.0;
  double lat = (y / 20037508.34) * 180.0;
  lat = 180.0 / M_PI * (2.0 * std::atan(std::exp(lat * M_PI / 180.0)) - M_PI / 2.0);
  return {lat, lon};
}

// 1x1 km square tile which can be visited for coverage
#define SQUADRAT_TILE_SIZE 1.0  // in km
// EPSG:3857 (Web Mercator) inflates east-west/north-south distances by 1/cos(lat) away from the
// equator, so a fixed-size square in projected meters isn't a real-world square everywhere.
// Scale the tile size so it measures exactly SQUADRAT_TILE_SIZE km at the reference latitude below.
#define SQUADRAT_REFERENCE_LATITUDE_DEG 55.6161  // Lejre, Denmark
static const double SQUADRAT_TILE_SIZE_METERS =
    SQUADRAT_TILE_SIZE * 1000.0 / std::cos(SQUADRAT_REFERENCE_LATITUDE_DEG * M_PI / 180.0);

struct SquadratTile {
  int squadrat_x{};  // X index of the tile in the grid EPSG 3857
  int squadrat_y{};  // Y index of the tile in the grid EPSG 3857

  bool operator<(const SquadratTile& other) const {
    return std::tie(squadrat_x, squadrat_y) < std::tie(other.squadrat_x, other.squadrat_y);
  }

  // Generate the tile from a point (lat, lon) in degrees
  SquadratTile(double lat, double lon) {
    // Convert to EPSG 3857 meters
    auto [x, y] = latLonToMeters(lat, lon);
    squadrat_x = static_cast<int>(std::floor(x / SQUADRAT_TILE_SIZE_METERS));
    squadrat_y = static_cast<int>(std::floor(y / SQUADRAT_TILE_SIZE_METERS));
  }

  // Get the bounding box of the tile in lat/lon degrees
  std::pair<Coordinate, Coordinate> getBBox() const {
    double min_x = squadrat_x * SQUADRAT_TILE_SIZE_METERS;
    double min_y = squadrat_y * SQUADRAT_TILE_SIZE_METERS;
    double max_x = (squadrat_x + 1) * SQUADRAT_TILE_SIZE_METERS;
    double max_y = (squadrat_y + 1) * SQUADRAT_TILE_SIZE_METERS;
    auto [min_lat, min_lon] = metersToLatLon(min_x, min_y);
    auto [max_lat, max_lon] = metersToLatLon(max_x, max_y);
    return {Coordinate{min_lat, min_lon}, Coordinate{max_lat, max_lon}};
  }
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SquadratTile, squadrat_x, squadrat_y)

class Tile {
 public:
  Tile(int z, int x, int y) : image_data_(256, 256, CV_8UC4, cv::Scalar(0, 0, 0, 0)), z_(z), x_(x), y_(y) {
    min_lon_ = x_ / std::pow(2.0, z_) * 360.0 - 180.0;
    max_lon_ = (x_ + 1) / std::pow(2.0, z_) * 360.0 - 180.0;
    min_lat_ = std::atan(std::sinh(M_PI * (1 - 2 * (y_ + 1) / std::pow(2.0, z_)))) * 180.0 / M_PI;
    max_lat_ = std::atan(std::sinh(M_PI * (1 - 2 * y_ / std::pow(2.0, z_)))) * 180.0 / M_PI;
    // fprintf(stderr, "Tile z=%d, x=%d, y=%d: min_lat=%.6f, max_lat=%.6f, min_lon=%.6f, max_lon=%.6f\n", z_, x_, y_,
    //         min_lat_, max_lat_, min_lon_, max_lon_);
  }

  void paint(const std::unordered_map<std::size_t, WaySegment>& way_segments) {
    int count = 0;
    for (const auto& [key, segment] : way_segments) {
      if (!segment.inBBox(min_lat_, min_lon_, max_lat_, max_lon_)) {
        continue;
      }
      count++;
      for (size_t i = 1; i < segment.geometry.size(); ++i) {
        const auto& p1 = segment.geometry[i - 1];
        const auto& p2 = segment.geometry[i];
        // fprintf(stderr, "Painting way %s (edge %llu) segment from (%.6f, %.6f) to (%.6f, %.6f)\n",
        //         segment.way_id.c_str(), static_cast<unsigned long long>(segment.edge_id), p1.lat, p1.lon, p2.lat,
        //         p2.lon);
        int x1 = static_cast<int>((p1.lon - min_lon_) / (max_lon_ - min_lon_) * 256);
        int y1 = static_cast<int>((max_lat_ - p1.lat) / (max_lat_ - min_lat_) * 256);
        int x2 = static_cast<int>((p2.lon - min_lon_) / (max_lon_ - min_lon_) * 256);
        int y2 = static_cast<int>((max_lat_ - p2.lat) / (max_lat_ - min_lat_) * 256);
        int alpha = std::min(255, segment.traversal_count * 10 + 50);  // Adjust the multiplier for desired opacity
        cv::line(image_data_, cv::Point(x1, y1), cv::Point(x2, y2), cv::Scalar(255, 0, 0, 255), 2, cv::LINE_AA);
      }
    }
    // fprintf(stderr, "Painted %d way segments on tile z=%d, x=%d, y=%d\n", count, z_, x_, y_);
  }

  // Draw the full squadrat grid (every 1km line), regardless of which tiles were visited
  void paintGrid() {
    if (z_ < 11) {
      return;  // Only draw grid for zoom levels 11 and above
    }
    int alpha = 100;  // Adjust the alpha value for desired opacity
    if (z_ >= 14) {
      alpha = 200;  // Make grid lines more visible at higher zoom levels
    }
    auto [min_x, min_y] = latLonToMeters(min_lat_, min_lon_);
    auto [max_x, max_y] = latLonToMeters(max_lat_, max_lon_);

    long long x_start = static_cast<long long>(std::floor(min_x / SQUADRAT_TILE_SIZE_METERS));
    long long x_end = static_cast<long long>(std::ceil(max_x / SQUADRAT_TILE_SIZE_METERS));
    cv::Scalar grey(20, 20, 20, alpha);
    for (long long i = x_start; i <= x_end; ++i) {
      auto [lat, lon] = metersToLatLon(i * SQUADRAT_TILE_SIZE_METERS, 0.0);
      int px = static_cast<int>((lon - min_lon_) / (max_lon_ - min_lon_) * 256);
      cv::line(image_data_, cv::Point(px, 0), cv::Point(px, 256), grey, 1, cv::LINE_AA);
    }

    long long y_start = static_cast<long long>(std::floor(min_y / SQUADRAT_TILE_SIZE_METERS));
    long long y_end = static_cast<long long>(std::ceil(max_y / SQUADRAT_TILE_SIZE_METERS));
    for (long long j = y_start; j <= y_end; ++j) {
      auto [lat, lon] = metersToLatLon(0.0, j * SQUADRAT_TILE_SIZE_METERS);
      int py = static_cast<int>((max_lat_ - lat) / (max_lat_ - min_lat_) * 256);
      cv::line(image_data_, cv::Point(0, py), cv::Point(256, py), grey, 1, cv::LINE_AA);
    }
  }

  void paint(const std::set<SquadratTile>& squadrat_tiles) {
    for (const auto& tile : squadrat_tiles) {
      auto [min_coord, max_coord] = tile.getBBox();
      int x1 = static_cast<int>((min_coord.lon - min_lon_) / (max_lon_ - min_lon_) * 256);
      int y1 = static_cast<int>((max_lat_ - max_coord.lat) / (max_lat_ - min_lat_) * 256);
      int x2 = static_cast<int>((max_coord.lon - min_lon_) / (max_lon_ - min_lon_) * 256);
      int y2 = static_cast<int>((max_lat_ - min_coord.lat) / (max_lat_ - min_lat_) * 256);
      cv::rectangle(image_data_, cv::Point(x1, y1), cv::Point(x2, y2), cv::Scalar(0, 255, 0, 100), cv::FILLED);
    }
  }

  cv::Mat getImage() const { return image_data_; }

  int getZ() const { return z_; }
  int getX() const { return x_; }
  int getY() const { return y_; }

  double minLon() const { return min_lon_; }
  double maxLon() const { return max_lon_; }
  double minLat() const { return min_lat_; }
  double maxLat() const { return max_lat_; }

 private:
  cv::Mat image_data_;
  int z_;
  int x_;
  int y_;
  double min_lat_;
  double max_lat_;
  double min_lon_;
  double max_lon_;
};

struct MatchedRoute {
  Route route;
  std::vector<WaySegment> way_segments;
  std::set<SquadratTile> squadrat_tiles;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MatchedRoute, route, way_segments, squadrat_tiles)

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
    request["trace_options"] = {{"search_radius", 50}, {"gps_accuracy", 20}};
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

    for (const auto& coord : route) {
      matched_route.squadrat_tiles.emplace(coord.lat, coord.lon);
    }

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
    matched_routes.reserve(routes.size());
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
      // if (matched_count >= 300) {
      //   break;
      // }
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

class TileKey {
 public:
  TileKey(int z, int x, int y) : z_(z), x_(x), y_(y) {}

  bool operator==(const TileKey& other) const { return z_ == other.z_ && x_ == other.x_ && y_ == other.y_; }

  struct Hash {
    std::size_t operator()(const TileKey& key) const {
      return std::hash<int>()(key.z_) ^ std::hash<int>()(key.x_) ^ std::hash<int>()(key.y_);
    }
  };

 private:
  int z_;
  int x_;
  int y_;
};

class TileGenerator {
 public:
  TileGenerator(const std::vector<MatchedRoute>& matched_routes) {
    for (const auto& matched_route : matched_routes) {
      squadrat_tiles_.insert(matched_route.squadrat_tiles.begin(), matched_route.squadrat_tiles.end());

      for (const auto& segment : matched_route.way_segments) {
        auto it = traversal_counts_.find(segment.edge_id);
        if (it == traversal_counts_.end()) {
          traversal_counts_[segment.edge_id] = segment;
        } else {
          it->second.traversal_count += 1;
        }
      }
    }
  }

  Tile generateTile(int z, int x, int y) {
    Tile tile(z, x, y);
    // tile.paint(traversal_counts_);
    tile.paintGrid();
    tile.paint(squadrat_tiles_);
    return tile;
  }

  Tile getTile(int z, int x, int y) {
    TileKey key(z, x, y);
    {
      std::lock_guard<std::mutex> lock(cache_mutex_);
      auto it = tile_cache_.find(key);
      if (it != tile_cache_.end()) {
        return it->second;
      }
    }

    // Generate the tile
    Tile tile = generateTile(z, x, y);

    // Cache the generated tile
    {
      std::lock_guard<std::mutex> lock(cache_mutex_);
      tile_cache_.emplace(key, tile);
    }

    return tile;
  }

 private:
  static std::mutex cache_mutex_;
  std::unordered_map<TileKey, Tile, TileKey::Hash> tile_cache_;
  std::unordered_map<std::size_t, WaySegment> traversal_counts_;
  std::set<SquadratTile> squadrat_tiles_;
};

std::mutex TileGenerator::cache_mutex_;

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

#if !RESTORE
  auto routes = load_all_routes();
  RouteMatcher matcher;
  auto matched_routes = matcher.matchAllRoutes(routes);
  nlohmann::json heatmap_json = matched_routes;
  std::ofstream heatmap_file(HEATMAP_FILE);
  if (!heatmap_file.is_open()) {
    throw std::runtime_error("Failed to open heatmap file for writing: " HEATMAP_FILE);
  }
  heatmap_file << heatmap_json.dump(2);
  heatmap_file.close();
#else
  std::ifstream heatmap_file(HEATMAP_FILE);
  if (!heatmap_file.is_open()) {
    throw std::runtime_error("Failed to open heatmap file for reading: " HEATMAP_FILE);
  }
  TimerLog restore_timer("Loading matched routes from heatmap.json");
  nlohmann::json heatmap_json;
  heatmap_file >> heatmap_json;
  std::vector<MatchedRoute> matched_routes = heatmap_json.get<std::vector<MatchedRoute> >();
  restore_timer.stop();
#endif

  TileGenerator tile_generator(matched_routes);

  // Main route planning endpoint
  svr.Get(R"(/tiles/(\d+)/(\d+)/(\d+).png)", [&tile_generator](const httplib::Request& req, httplib::Response& res) {
    // Heatmap XYZ tile request from url like /tiles/{z}/{x}/{y}.png
    for (const auto& param : req.path_params) {
      fprintf(stderr, "Path param: %s = %s\n", param.first.c_str(), param.second.c_str());
    }
    int z = std::stoi(req.matches[1]);
    int x = std::stoi(req.matches[2]);
    int y = std::stoi(req.matches[3]);
    // fprintf(stderr, "Received tile request for z=%d, x=%d, y=%d\n", z, x, y);

    // For now just return a placeholder PNG image (256x256 pixel)
    Tile tile = tile_generator.getTile(z, x, y);
    cv::Mat png_data = tile.getImage();
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
