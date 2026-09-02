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
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <opencv2/core/mat.hpp>
#include <opencv4/opencv2/opencv.hpp>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "httplib.h"
#include "json.hpp"

#define VALHALLA_CONFIG_FILE "/home/christian/git/witx-heatmap/data/routing/valhalla_data/valhalla.json"
#define RESTORE 1

// Christian
#define ROUTE_FILE "/media/christian/Data/Backup/strava/strava_christian_simple.geojson"
#define HEATMAP_FILE "/home/christian/git/witx-heatmap/data/heatmap.json"

// Thomas
// #define ROUTE_FILE "/media/christian/Data/Backup/strava/strava_thomas_simple.geojson"
// #define HEATMAP_FILE "/home/christian/git/witx-heatmap/data/heatmap_thomas.json"

struct Coordinate {
  double lat{};
  double lon{};
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Coordinate, lat, lon)

// using Route = std::vector<Coordinate>;

struct Route {
  std::string link;
  std::string name;
  std::vector<Coordinate> route;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Route, link, name, route)

struct WaySegment {
  std::string way_id;
  uint64_t edge_id{};     // Valhalla's internal directed-edge id: stable identity for the exact
                          // maximal way-segment traversed, safe to use as a traversal-count key
  Route geometry;         // The portion of the way's geometry actually traversed by the route
  int traversal_count{};  // How many times this way-segment was traversed across all matched routes

  bool inBBox(double min_lat, double min_lon, double max_lat, double max_lon) const {
    for (int i = 0; i < geometry.route.size(); ++i) {
      const auto& coord = geometry.route[i];
      if (coord.lat >= min_lat && coord.lat <= max_lat && coord.lon >= min_lon && coord.lon <= max_lon) {
        return true;
      }
      if (i > 0) {
        const auto& prev_coord = geometry.route[i - 1];
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
// EPSG:3857 (Web Mercator) inflates east-west/north-south distances by 1/cos(lat) away from the
// equator, so a fixed-size square in projected meters isn't a real-world square everywhere.
// Scale the tile size so it measures exactly SQUADRAT_TILE_SIZE km at the reference latitude below.
#define SQUADRAT_REFERENCE_LATITUDE_DEG 55.6161  // Lejre, Denmark
double meter2size(double size) {
  return size / std::cos(SQUADRAT_REFERENCE_LATITUDE_DEG * M_PI / 180.0);
};

struct SquadratTile {
  int squadrat_x{};    // X index of the tile in the grid EPSG 3857
  int squadrat_y{};    // Y index of the tile in the grid EPSG 3857
  double tile_size{};  // Size of the tile in corrected meters

  bool operator<(const SquadratTile& other) const {
    return std::tie(squadrat_x, squadrat_y) < std::tie(other.squadrat_x, other.squadrat_y);
  }

  SquadratTile() = default;  // needed for JSON deserialization

  // Generate the tile from a point (lat, lon) in degrees
  SquadratTile(double lat, double lon, double tile_size) {
    // Convert to EPSG 3857 meters
    auto [x, y] = latLonToMeters(lat, lon);
    squadrat_x = static_cast<int>(std::floor(x / tile_size));
    squadrat_y = static_cast<int>(std::floor(y / tile_size));
    this->tile_size = tile_size;
  }

  // Get the bounding box of the tile in lat/lon degrees
  std::pair<Coordinate, Coordinate> getBBox() const {
    double min_x = squadrat_x * tile_size;
    double min_y = squadrat_y * tile_size;
    double max_x = (squadrat_x + 1) * tile_size;
    double max_y = (squadrat_y + 1) * tile_size;
    auto [min_lat, min_lon] = metersToLatLon(min_x, min_y);
    auto [max_lat, max_lon] = metersToLatLon(max_x, max_y);
    return {Coordinate{min_lat, min_lon}, Coordinate{max_lat, max_lon}};
  }
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SquadratTile, squadrat_x, squadrat_y, tile_size)

struct MatchedRoute {
  Route route;
  std::vector<WaySegment> way_segments;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MatchedRoute, route, way_segments)

static double haversineDistance(const Coordinate& coord1, const Coordinate& coord2) {
  constexpr double EARTH_RADIUS_KM = 6371.0;
  double lat1_rad = coord1.lat * M_PI / 180.0;
  double lon1_rad = coord1.lon * M_PI / 180.0;
  double lat2_rad = coord2.lat * M_PI / 180.0;
  double lon2_rad = coord2.lon * M_PI / 180.0;

  double dlat = lat2_rad - lat1_rad;
  double dlon = lon2_rad - lon1_rad;

  double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
             std::cos(lat1_rad) * std::cos(lat2_rad) * std::sin(dlon / 2) * std::sin(dlon / 2);
  double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

  return EARTH_RADIUS_KM * c;
}

class AlphaShape {
 public:
  // alpha is the max circumradius (in meters) a Delaunay triangle may have to stay in the shape;
  // smaller alpha follows the point cloud's concavities more tightly, larger alpha tends to the convex hull.
  AlphaShape(const std::vector<MatchedRoute>& matched_routes, double alpha) : alpha_(alpha) {
    constexpr Coordinate START_COORD{55.59784, 11.97298};
    constexpr double TOLERANCE_KM = 0.2;
    for (const auto& matched_route : matched_routes) {
      if (matched_route.route.empty()) {
        continue;
      }
      const auto& start_coord = matched_route.route.front();
      double distance_km = haversineDistance(start_coord, START_COORD);
      if (distance_km > TOLERANCE_KM) {
        continue;
      }
      const auto& end_coord = matched_route.route.back();
      double end_distance_km = haversineDistance(end_coord, START_COORD);
      if (end_distance_km > TOLERANCE_KM) {
        continue;
      }

      matched_routes_.push_back(matched_route);

      for (const auto& coord : matched_route.route.route) {
        points_.emplace_back(coord.lat, coord.lon);
      }
    }
    compute();
  }

  const std::vector<std::pair<Coordinate, Coordinate>>& getBoundaryEdges() const {
    return boundary_edges_;
  }

  void compute() {
    boundary_edges_.clear();
    if (points_.size() < 3) {
      return;
    }

    // Project to planar meters so the alpha radius test means the same thing everywhere.
    std::vector<cv::Point2f> projected;
    projected.reserve(points_.size());
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    for (const auto& [lat, lon] : points_) {
      auto [x, y] = latLonToMeters(lat, lon);
      min_x = std::min(min_x, x);
      min_y = std::min(min_y, y);
      max_x = std::max(max_x, x);
      max_y = std::max(max_y, y);
      projected.emplace_back(static_cast<float>(x), static_cast<float>(y));
    }

    // Subdiv2D requires a rect that strictly contains every inserted point.
    cv::Rect2f bounds(static_cast<float>(min_x - 1.0), static_cast<float>(min_y - 1.0),
                      static_cast<float>(max_x - min_x + 2.0), static_cast<float>(max_y - min_y + 2.0));
    cv::Subdiv2D subdiv(bounds);
    for (const auto& p : projected) {
      subdiv.insert(p);
    }

    std::vector<cv::Vec6f> triangles;
    subdiv.getTriangleList(triangles);

    // An edge that borders exactly one alpha-valid triangle lies on the boundary of the shape.
    auto edge_key = [](cv::Point2f a, cv::Point2f b) {
      if (std::make_pair(a.x, a.y) > std::make_pair(b.x, b.y)) {
        std::swap(a, b);
      }
      return std::make_tuple(a.x, a.y, b.x, b.y);
    };
    std::map<std::tuple<float, float, float, float>, int> edge_counts;
    std::map<std::tuple<float, float, float, float>, std::pair<cv::Point2f, cv::Point2f>> edge_points;

    for (const auto& t : triangles) {
      cv::Point2f p1(t[0], t[1]);
      cv::Point2f p2(t[2], t[3]);
      cv::Point2f p3(t[4], t[5]);

      // Subdiv2D generates extra triangles connecting to its bounding rect; ignore those.
      if (!bounds.contains(p1) || !bounds.contains(p2) || !bounds.contains(p3)) {
        continue;
      }

      if (circumradius(p1, p2, p3) > alpha_) {
        continue;
      }

      for (const auto& [a, b] : {std::make_pair(p1, p2), std::make_pair(p2, p3), std::make_pair(p3, p1)}) {
        auto key = edge_key(a, b);
        edge_counts[key]++;
        edge_points[key] = {a, b};
      }
    }

    for (const auto& [key, count] : edge_counts) {
      if (count != 1) {
        continue;
      }
      const auto& [a, b] = edge_points[key];
      auto [lat1, lon1] = metersToLatLon(a.x, a.y);
      auto [lat2, lon2] = metersToLatLon(b.x, b.y);
      boundary_edges_.emplace_back(Coordinate{lat1, lon1}, Coordinate{lat2, lon2});
    }
  }

  const std::vector<MatchedRoute>& getMatchedRoutes() const {
    return matched_routes_;
  }

 private:
  // Circumradius of triangle (a, b, c); returns +inf for degenerate (near-zero-area) triangles.
  static double circumradius(cv::Point2f a, cv::Point2f b, cv::Point2f c) {
    double ab = cv::norm(a - b);
    double bc = cv::norm(b - c);
    double ca = cv::norm(c - a);
    double s = (ab + bc + ca) / 2.0;
    double area = std::sqrt(std::max(0.0, s * (s - ab) * (s - bc) * (s - ca)));
    if (area < 1e-9) {
      return std::numeric_limits<double>::infinity();
    }
    return (ab * bc * ca) / (4.0 * area);
  }

  double alpha_;
  std::vector<std::pair<double, double>> points_;
  std::vector<std::pair<Coordinate, Coordinate>> boundary_edges_;
  std::vector<MatchedRoute> matched_routes_;
};

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

  void paintMatches(const std::unordered_map<std::size_t, WaySegment>& way_segments) {
    int count = 0;
    for (const auto& [key, segment] : way_segments) {
      if (!segment.inBBox(min_lat_, min_lon_, max_lat_, max_lon_)) {
        continue;
      }
      count++;
      for (size_t i = 1; i < segment.geometry.route.size(); ++i) {
        const auto& p1 = segment.geometry.route[i - 1];
        const auto& p2 = segment.geometry.route[i];
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

  void paintRoute(const std::vector<MatchedRoute>& matched_routes) {
    int count = 0;
    for (const auto& matched_route : matched_routes) {
      for (size_t i = 1; i < matched_route.route.route.size(); ++i) {
        const auto& p1 = matched_route.route.route[i - 1];
        const auto& p2 = matched_route.route.route[i];
        // fprintf(stderr, "Painting way %s (edge %llu) segment from (%.6f, %.6f) to (%.6f, %.6f)\n",
        //         segment.way_id.c_str(), static_cast<unsigned long long>(segment.edge_id), p1.lat, p1.lon, p2.lat,
        //         p2.lon);
        int x1 = static_cast<int>((p1.lon - min_lon_) / (max_lon_ - min_lon_) * 256);
        int y1 = static_cast<int>((max_lat_ - p1.lat) / (max_lat_ - min_lat_) * 256);
        int x2 = static_cast<int>((p2.lon - min_lon_) / (max_lon_ - min_lon_) * 256);
        int y2 = static_cast<int>((max_lat_ - p2.lat) / (max_lat_ - min_lat_) * 256);
        int alpha = 255;  // Adjust the multiplier for desired opacity
        cv::line(image_data_, cv::Point(x1, y1), cv::Point(x2, y2), cv::Scalar(255, 0, 0, 255), 2, cv::LINE_AA);
      }
    }
    // fprintf(stderr, "Painted %d way segments on tile z=%d, x=%d, y=%d\n", count, z_, x_, y_);
  }

  // Draw the full squadrat grid (every 1km line), regardless of which tiles were visited
  void paintGrid(double tile_size) {
    if (z_ < 11) {
      return;  // Only draw grid for zoom levels 11 and above
    }
    int alpha = 100;  // Adjust the alpha value for desired opacity
    if (z_ >= 14) {
      alpha = 200;  // Make grid lines more visible at higher zoom levels
    }
    auto [min_x, min_y] = latLonToMeters(min_lat_, min_lon_);
    auto [max_x, max_y] = latLonToMeters(max_lat_, max_lon_);

    long long x_start = static_cast<long long>(std::floor(min_x / tile_size));
    long long x_end = static_cast<long long>(std::ceil(max_x / tile_size));
    cv::Scalar grey(20, 20, 20, alpha);
    for (long long i = x_start; i <= x_end; ++i) {
      auto [lat, lon] = metersToLatLon(i * tile_size, 0.0);
      int px = static_cast<int>((lon - min_lon_) / (max_lon_ - min_lon_) * 256);
      cv::line(image_data_, cv::Point(px, 0), cv::Point(px, 256), grey, 1, cv::LINE_AA);
    }

    long long y_start = static_cast<long long>(std::floor(min_y / tile_size));
    long long y_end = static_cast<long long>(std::ceil(max_y / tile_size));
    for (long long j = y_start; j <= y_end; ++j) {
      auto [lat, lon] = metersToLatLon(0.0, j * tile_size);
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

  void paint(const AlphaShape& alpha_shape) {
    paintRoute(alpha_shape.getMatchedRoutes());
    for (const auto& [p1, p2] : alpha_shape.getBoundaryEdges()) {
      int x1 = static_cast<int>((p1.lon - min_lon_) / (max_lon_ - min_lon_) * 256);
      int y1 = static_cast<int>((max_lat_ - p1.lat) / (max_lat_ - min_lat_) * 256);
      int x2 = static_cast<int>((p2.lon - min_lon_) / (max_lon_ - min_lon_) * 256);
      int y2 = static_cast<int>((max_lat_ - p2.lat) / (max_lat_ - min_lat_) * 256);
      cv::line(image_data_, cv::Point(x1, y1), cv::Point(x2, y2), cv::Scalar(255, 0, 255, 200), 3, cv::LINE_AA);
    }
  }

  cv::Mat getImage() const {
    return image_data_;
  }

  int getZ() const {
    return z_;
  }
  int getX() const {
    return x_;
  }
  int getY() const {
    return y_;
  }

  double minLon() const {
    return min_lon_;
  }
  double maxLon() const {
    return max_lon_;
  }
  double minLat() const {
    return min_lat_;
  }
  double maxLat() const {
    return max_lat_;
  }

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
    route.link = feature["properties"]["link"];
    route.name = feature["properties"]["name"];
    for (const auto& point : geometry["coordinates"]) {
      // GeoJSON coordinates are [lon, lat]
      route.route.push_back(Coordinate{point[1].get<double>(), point[0].get<double>()});
    }
    routes.push_back(std::move(route));
  }

  return routes;
}
// Using chrono timer
class TimerLog {
 public:
  TimerLog(const std::string& name) : name_(name), start_(std::chrono::high_resolution_clock::now()) {
  }

  ~TimerLog() {
    stop();
  }

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
  MatchedRoute matchRoute(const Route& route) {
    nlohmann::json request;
    request["costing"] = "pedestrian";
    request["shape_match"] = "map_snap";
    request["trace_options"] = {{"search_radius", 50}, {"gps_accuracy", 20}};
    for (const auto& coord : route.route) {
      request["shape"].push_back({{"lat", coord.lat}, {"lon", coord.lon}});
    }

    MatchedRoute matched_route;
    matched_route.route = route;

    valhalla::Api api;
    std::string json_str;
    try {
      json_str = actor_->trace_route(request.dump(), nullptr, &api);
    } catch (const std::exception& e) {
      fprintf(stderr, "Valhalla trace_route failed: %s\n", e.what());
      return matched_route;
    }

    // print the raw json result to stderr
    // fprintf(stderr, "Valhalla trace_route result: %s\n", json_str.c_str());

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
          Route geometry;
          geometry.route.insert(geometry.route.end(), leg_shape.begin() + begin, leg_shape.begin() + end + 1);
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
      if (matched_route.way_segments.empty()) {
        fprintf(stderr, "Route matching failed for: %s, %s\n", route.name.c_str(), route.link.c_str());
      }
      matched_routes.push_back(std::move(matched_route));
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
  TileKey(int z, int x, int y) : z_(z), x_(x), y_(y) {
  }

  bool operator==(const TileKey& other) const {
    return z_ == other.z_ && x_ == other.x_ && y_ == other.y_;
  }

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
  virtual Tile generateTile(int z, int x, int y) = 0;

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
  std::mutex cache_mutex_;
  std::unordered_map<TileKey, Tile, TileKey::Hash> tile_cache_;
};

class AlphaShapeTileGenerator : public TileGenerator {
 public:
  AlphaShapeTileGenerator(const std::vector<MatchedRoute>& matched_routes, double alpha)
      : alpha_shape_(matched_routes, alpha) {
  }

  Tile generateTile(int z, int x, int y) override {
    Tile tile(z, x, y);
    tile.paint(alpha_shape_);
    return tile;
  }

 private:
  AlphaShape alpha_shape_;
};

class SquadratTileGenerator : public TileGenerator {
 public:
  SquadratTileGenerator(const std::vector<MatchedRoute>& matched_routes, double tile_size_raw)
      : tile_size_(meter2size(tile_size_raw)) {
    for (const auto& matched_route : matched_routes) {
      for (const auto& coordinate : matched_route.route.route) {
        squadrat_tiles_.emplace(coordinate.lat, coordinate.lon, tile_size_);
      }
    }
  }

  Tile generateTile(int z, int x, int y) override {
    Tile tile(z, x, y);
    tile.paint(squadrat_tiles_);
    tile.paintGrid(tile_size_);
    return tile;
  }

 private:
  std::set<SquadratTile> squadrat_tiles_;
  double tile_size_;
};

class TraversalTileGenerator : public TileGenerator {
 public:
  TraversalTileGenerator(const std::vector<MatchedRoute>& matched_routes) {
    for (const auto& matched_route : matched_routes) {
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

  Tile generateTile(int z, int x, int y) override {
    Tile tile(z, x, y);
    tile.paintMatches(traversal_counts_);
    return tile;
  }

 private:
  std::unordered_map<std::size_t, WaySegment> traversal_counts_;
};

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
  std::vector<MatchedRoute> matched_routes = heatmap_json.get<std::vector<MatchedRoute>>();
  restore_timer.stop();
#endif

  std::unordered_map<int, std::unique_ptr<AlphaShapeTileGenerator>> alpha_shapes;
  {
    TimerLog alpha_shapes_timer("Generating alpha shapes for radius 4000");
    alpha_shapes.emplace(4000, std::make_unique<AlphaShapeTileGenerator>(matched_routes, static_cast<double>(4000)));
  }
  {
    TimerLog alpha_shapes_timer("Generating alpha shapes for radius 7000");
    alpha_shapes.emplace(7000, std::make_unique<AlphaShapeTileGenerator>(matched_routes, static_cast<double>(7000)));
  }
  {
    TimerLog alpha_shapes_timer("Generating alpha shapes for radius 10000");
    alpha_shapes.emplace(10000, std::make_unique<AlphaShapeTileGenerator>(matched_routes, static_cast<double>(10000)));
  }
  std::mutex alpha_shapes_mutex;

  // Main route planning endpoint
  svr.Get(R"(/coverage/(\d+)/(\d+)/(\d+).png)", [&matched_routes, &alpha_shapes, &alpha_shapes_mutex](
                                                    const httplib::Request& req, httplib::Response& res) {
    // Heatmap XYZ tile request from url like /tiles/{z}/{x}/{y}.png
    for (const auto& param : req.path_params) {
      fprintf(stderr, "Path param: %s = %s\n", param.first.c_str(), param.second.c_str());
    }
    int z = std::stoi(req.matches[1]);
    int x = std::stoi(req.matches[2]);
    int y = std::stoi(req.matches[3]);
    std::string user = req.has_param("user") ? req.get_param_value("user") : "";
    int radius = req.has_param("radius") ? std::stoi(req.get_param_value("radius")) : 0;
    // fprintf(stderr, "Received tile request for z=%d, x=%d, y=%d\n", z, x, y);

    // For now just return a placeholder PNG image (256x256 pixel)
    auto it = alpha_shapes.find(radius);
    if (it == alpha_shapes.end()) {
      std::lock_guard<std::mutex> lock(alpha_shapes_mutex);
      it = alpha_shapes
               .emplace(radius, std::make_unique<AlphaShapeTileGenerator>(matched_routes, static_cast<double>(radius)))
               .first;
    }
    Tile tile = it->second->getTile(z, x, y);
    cv::Mat png_data = tile.getImage();
    std::vector<unsigned char> buffer;
    cv::imencode(".png", png_data, buffer);
    res.set_content(reinterpret_cast<const char*>(buffer.data()), buffer.size(), "image/png");
  });

  std::unordered_map<int, std::unique_ptr<SquadratTileGenerator>> squadrat_tile_generators;
  {
    TimerLog squadrat_tiles_timer("Generating squadrats for radius 4000");
    squadrat_tile_generators.emplace(
        4000, std::make_unique<SquadratTileGenerator>(matched_routes, static_cast<double>(1000)));
  }
  {
    TimerLog squadrat_tiles_timer("Generating squadrats for radius 7000");
    squadrat_tile_generators.emplace(
        1600, std::make_unique<SquadratTileGenerator>(matched_routes, static_cast<double>(1600)));
  }
  {
    TimerLog squadrat_tiles_timer("Generating squadrats for radius 10000");
    squadrat_tile_generators.emplace(
        10000, std::make_unique<SquadratTileGenerator>(matched_routes, static_cast<double>(5000)));
  }
  std::mutex squadrat_tiles_mutex;

  svr.Get(R"(/squadrat/(\d+)/(\d+)/(\d+).png)", [&matched_routes, &squadrat_tile_generators, &squadrat_tiles_mutex](
                                                    const httplib::Request& req, httplib::Response& res) {
    // Heatmap XYZ tile request from url like /tiles/{z}/{x}/{y}.png
    for (const auto& param : req.path_params) {
      fprintf(stderr, "Path param: %s = %s\n", param.first.c_str(), param.second.c_str());
    }
    int z = std::stoi(req.matches[1]);
    int x = std::stoi(req.matches[2]);
    int y = std::stoi(req.matches[3]);
    std::string user = req.has_param("user") ? req.get_param_value("user") : "";
    int radius = req.has_param("radius") ? std::stoi(req.get_param_value("radius")) : 0;
    // fprintf(stderr, "Received tile request for z=%d, x=%d, y=%d\n", z, x, y);

    // For now just return a placeholder PNG image (256x256 pixel)
    auto it = squadrat_tile_generators.find(radius);
    if (it == squadrat_tile_generators.end()) {
      std::lock_guard<std::mutex> lock(squadrat_tiles_mutex);
      it = squadrat_tile_generators
               .emplace(radius, std::make_unique<SquadratTileGenerator>(matched_routes, static_cast<double>(radius)))
               .first;
    }
    Tile tile = it->second->getTile(z, x, y);
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
