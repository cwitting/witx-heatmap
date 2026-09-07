// #include <omp.h>

#include <valhalla/config.h>
#include <valhalla/tyr/actor.h>

#include <algorithm>
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
#include <shared_mutex>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "heatmap.pb.h"
#include "httplib.h"
#include "json.hpp"

#define VALHALLA_CONFIG_FILE "/home/christian/git/witx-heatmap/data/routing/valhalla_data/valhalla.json"
#define RESTORE 0

const std::set<std::string> activity_id_blacklist = {
    "5055464955",
};

// Christian
#define ROUTE_FILE "/media/christian/Data/Backup/strava/strava_christian_full.geojson"
#define HEATMAP_FILE "/media/christian/Data/Backup/strava/heatmap_christian.pb"
#define DEFAULT_INGEST_FOLDER "/home/christian/git/witx-heatmap/data/ingest_christian"

// Thomas
// #define ROUTE_FILE "/media/christian/Data/Backup/strava/strava_thomas_full.geojson"
// #define HEATMAP_FILE "/media/christian/Data/Backup/strava/heatmap_thomas.pb"
// #define DEFAULT_INGEST_FOLDER "/home/christian/git/witx-heatmap/data/ingest_thomas"

// Nikolaj
// #define ROUTE_FILE "/media/christian/Data/Backup/strava/strava_nikolaj_full.geojson"
// #define HEATMAP_FILE "/home/christian/git/witx-heatmap/data/heatmap_nikolaj.pb"
// #define DEFAULT_INGEST_FOLDER "/home/christian/git/witx-heatmap/data/ingest_nikolaj"

struct Coordinate {
  double lat{};
  double lon{};
  bool isValid() const {
    return lat >= -90.0 && lat <= 90.0 && lon >= -180.0 && lon <= 180.0;
  }
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Coordinate, lat, lon)

// Simplify a coordinate sequence by dropping points too close to the last kept point (lon/lat distance only).
static std::vector<Coordinate> simplifyCoordinates(const std::vector<Coordinate>& coords, double tolerance = 0.001) {
  if (coords.size() < 2) {
    return coords;
  }

  std::vector<Coordinate> simplified{coords[0]};
  for (std::size_t i = 1; i < coords.size(); ++i) {
    const auto& last = simplified.back();
    const auto& coord = coords[i];
    double dist =
        std::sqrt((coord.lon - last.lon) * (coord.lon - last.lon) + (coord.lat - last.lat) * (coord.lat - last.lat));
    if (dist > tolerance) {
      simplified.push_back(coord);
    }
  }

  return simplified;
}

constexpr double MILLISECONDS_PER_DAY = 24 * 60 * 60 * 1000;

static double millisecondsNow() {
  return std::chrono::duration<double, std::milli>(std::chrono::system_clock::now().time_since_epoch()).count();
}

static std::vector<cv::Scalar> default_colors = {
    cv::Scalar(0, 255, 0, 130),    // green
    cv::Scalar(0, 255, 255, 130),  // yellow
    cv::Scalar(0, 0, 255, 130)     // red
};

cv::Scalar color_map(double value, std::vector<cv::Scalar> colors) {
  if (colors.empty()) {
    return cv::Scalar(0, 0, 0, 200);  // default to black if no colors provided
  }
  value = std::clamp(value, 0.0, 1.0);
  std::size_t lower_index = static_cast<std::size_t>(value * (colors.size() - 1));
  std::size_t upper_index = std::min(lower_index + 1, colors.size() - 1);
  double t = value * (colors.size() - 1) - lower_index;
  cv::Scalar lower_color = colors[lower_index];
  cv::Scalar upper_color = colors[upper_index];
  return lower_color * (1.0 - t) + upper_color * t;
}

struct Route {
  std::string link;
  std::string activity_id;
  std::string name;
  std::string date;
  int date_format{};
  std::vector<Coordinate> route;
  double getMilliseconds() const {
    // Format is Sep 1, 2026, 12:52:03 PM
    std::tm tm{};
    if (date_format == 0) {
      if (strptime(date.c_str(), "%b %d, %Y, %I:%M:%S %p", &tm) == nullptr) {
        throw std::runtime_error("Failed to parse date: " + date);
      }
    } else {
      // Format: 2026-08-29T11:21:42Z
      std::istringstream ss(date);
      ss >> std::get_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
      if (ss.fail()) {
        throw std::runtime_error("Failed to parse date string: " + date);
      }
    }
    return static_cast<double>(std::mktime(&tm)) * 1000.0;
  }
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Route, link, activity_id, name, date, route)

struct WaySegment {
  std::string way_id;
  uint64_t edge_id{};     // Valhalla's internal directed-edge id: stable identity for the exact
                          // maximal way-segment traversed, safe to use as a traversal-count key
  Route geometry;         // The portion of the way's geometry actually traversed by the route
  int traversal_count{};  // How many times this way-segment was traversed across all matched routes
  double first_traversal_time{};
  double last_traversal_time{};

  void visit(double time) {
    if (first_traversal_time == 0 || time < first_traversal_time) {
      first_traversal_time = time;
    }
    if (time > last_traversal_time) {
      last_traversal_time = time;
    }
  }

  double getFirstVisitTime() const {
    return first_traversal_time;
  }

  double getLastVisitTime() const {
    return last_traversal_time;
  }

  double getFirstVisitAge() const {
    return std::max(0.0, millisecondsNow() - first_traversal_time) / MILLISECONDS_PER_DAY;
  }

  double getLastVisitAge() const {
    return std::max(0.0, millisecondsNow() - last_traversal_time) / MILLISECONDS_PER_DAY;
  }

  static constexpr double AGE_THRESHOLD = 365;

  cv::Scalar getFirstVisitColor() const {
    double age = std::min(1.0, getFirstVisitAge() / AGE_THRESHOLD);
    cv::Scalar color = color_map(age, default_colors);
    color[3] = 255;  // Ensure the alpha channel is set to 255
    return color;
  }

  cv::Scalar getLastVisitColor() const {
    double age = std::min(1.0, getLastVisitAge() / AGE_THRESHOLD);
    cv::Scalar color = color_map(age, default_colors);
    color[3] = 255;  // Ensure the alpha channel is set to 255
    return color;
  }

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

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(WaySegment, way_id, edge_id, geometry, traversal_count, first_traversal_time,
                                   last_traversal_time)

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

constexpr double HEX_SQRT3 = 1.7320508075688772;

// 1x1 km hexagonal tile which can be visited for coverage. Uses "pointy-top" axial coordinates
// (q, r), see https://www.redblobgames.com/grids/hexagons/ for the underlying math.
struct SquadratTile {
  int squadrat_q{};    // Axial q index of the hexagon in the grid EPSG 3857
  int squadrat_r{};    // Axial r index of the hexagon in the grid EPSG 3857
  double tile_size{};  // Circumradius (center to corner) of the hexagon in corrected meters
  mutable double first_visit_time{};
  mutable double last_visit_time{};

  bool operator<(const SquadratTile& other) const {
    return std::tie(squadrat_q, squadrat_r) < std::tie(other.squadrat_q, other.squadrat_r);
  }

  SquadratTile() = default;  // needed for JSON deserialization

  // Generate the tile from a point (lat, lon) in degrees
  SquadratTile(double lat, double lon, double tile_size) {
    // Convert to EPSG 3857 meters
    auto [x, y] = latLonToMeters(lat, lon);
    double qf = (HEX_SQRT3 / 3.0 * x - y / 3.0) / tile_size;
    double rf = (2.0 / 3.0 * y) / tile_size;
    cubeRound(qf, rf, squadrat_q, squadrat_r);
    this->tile_size = tile_size;
  }

  void visit(double time) const {
    if (first_visit_time == 0 || time < first_visit_time) {
      first_visit_time = time;
    }
    if (time > last_visit_time) {
      last_visit_time = time;
    }
  }

  double getFirstVisitTime() const {
    return first_visit_time;
  }

  double getLastVisitTime() const {
    return last_visit_time;
  }

  double getFirstVisitAge() const {
    return std::max(0.0, millisecondsNow() - first_visit_time) / MILLISECONDS_PER_DAY;
  }

  double getLastVisitAge() const {
    return std::max(0.0, millisecondsNow() - last_visit_time) / MILLISECONDS_PER_DAY;
  }

  static constexpr double AGE_THRESHOLD = 365;

  cv::Scalar getFirstVisitColor() const {
    double age = std::min(1.0, getFirstVisitAge() / AGE_THRESHOLD);
    // Age 0 = fresh (green), Age 1 = old (red)
    return color_map(age, default_colors);
  }

  cv::Scalar getLastVisitColor() const {
    double age = std::min(1.0, getLastVisitAge() / AGE_THRESHOLD);
    return color_map(age, default_colors);
  }

  // Get the 6 corners of the hexagon in lat/lon degrees
  std::vector<Coordinate> getVertices() const {
    return hexVertices(squadrat_q, squadrat_r, tile_size);
  }

  // Axial (q, r) hex center, in the same corrected EPSG 3857 meters as tile_size
  static std::pair<double, double> hexCenter(int q, int r, double tile_size) {
    double x = tile_size * (HEX_SQRT3 * q + HEX_SQRT3 / 2.0 * r);
    double y = tile_size * 1.5 * r;
    return {x, y};
  }

  static std::vector<Coordinate> hexVertices(int q, int r, double tile_size) {
    auto [cx, cy] = hexCenter(q, r, tile_size);
    std::vector<Coordinate> vertices;
    vertices.reserve(6);
    for (int i = 0; i < 6; ++i) {
      double angle = M_PI / 180.0 * (60.0 * i - 30.0);
      auto [lat, lon] = metersToLatLon(cx + tile_size * std::cos(angle), cy + tile_size * std::sin(angle));
      vertices.push_back(Coordinate{lat, lon});
    }
    return vertices;
  }

 private:
  // Round fractional cube coordinates (qf, rf, -qf-rf) to the nearest hex, fixing up whichever
  // component has the largest rounding error so q+r+s stays exactly zero.
  static void cubeRound(double qf, double rf, int& q, int& r) {
    double sf = -qf - rf;
    double q_round = std::round(qf);
    double r_round = std::round(rf);
    double s_round = std::round(sf);
    double q_diff = std::abs(q_round - qf);
    double r_diff = std::abs(r_round - rf);
    double s_diff = std::abs(s_round - sf);
    if (q_diff > r_diff && q_diff > s_diff) {
      q_round = -r_round - s_round;
    } else if (r_diff > s_diff) {
      r_round = -q_round - s_round;
    }
    q = static_cast<int>(q_round);
    r = static_cast<int>(r_round);
  }
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SquadratTile, squadrat_q, squadrat_r, tile_size, first_visit_time, last_visit_time)

struct MatchedRoute {
  Route route;
  std::vector<WaySegment> way_segments;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MatchedRoute, route, way_segments)

static witxheatmap::PRoute routeToProto(const Route& route) {
  witxheatmap::PRoute proto;
  proto.set_link(route.link);
  proto.set_activity_id(route.activity_id);
  proto.set_name(route.name);
  proto.set_date(route.date);
  for (const auto& coord : route.route) {
    auto* proto_coord = proto.add_route();
    proto_coord->set_lat(coord.lat);
    proto_coord->set_lon(coord.lon);
  }
  return proto;
}

static Route routeFromProto(const witxheatmap::PRoute& proto) {
  Route route;
  route.link = proto.link();
  route.activity_id = proto.activity_id();
  route.name = proto.name();
  route.date = proto.date();
  route.route.reserve(proto.route_size());
  for (const auto& proto_coord : proto.route()) {
    route.route.push_back(Coordinate{proto_coord.lat(), proto_coord.lon()});
  }
  return route;
}

static witxheatmap::PWaySegment waySegmentToProto(const WaySegment& segment) {
  witxheatmap::PWaySegment proto;
  proto.set_way_id(segment.way_id);
  proto.set_edge_id(segment.edge_id);
  *proto.mutable_geometry() = routeToProto(segment.geometry);
  proto.set_traversal_count(segment.traversal_count);
  proto.set_first_traversal_time(segment.first_traversal_time);
  proto.set_last_traversal_time(segment.last_traversal_time);
  return proto;
}

static WaySegment waySegmentFromProto(const witxheatmap::PWaySegment& proto) {
  WaySegment segment;
  segment.way_id = proto.way_id();
  segment.edge_id = proto.edge_id();
  segment.geometry = routeFromProto(proto.geometry());
  segment.traversal_count = proto.traversal_count();
  segment.first_traversal_time = proto.first_traversal_time();
  segment.last_traversal_time = proto.last_traversal_time();
  return segment;
}

static witxheatmap::PMatchedRoute matchedRouteToProto(const MatchedRoute& matched_route) {
  witxheatmap::PMatchedRoute proto;
  *proto.mutable_route() = routeToProto(matched_route.route);
  for (const auto& segment : matched_route.way_segments) {
    *proto.add_way_segments() = waySegmentToProto(segment);
  }
  return proto;
}

static MatchedRoute matchedRouteFromProto(const witxheatmap::PMatchedRoute& proto) {
  MatchedRoute matched_route;
  matched_route.route = routeFromProto(proto.route());
  matched_route.way_segments.reserve(proto.way_segments_size());
  for (const auto& proto_segment : proto.way_segments()) {
    matched_route.way_segments.push_back(waySegmentFromProto(proto_segment));
  }
  return matched_route;
}

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
  AlphaShape() = default;
  // alpha is the max circumradius (in meters) a Delaunay triangle may have to stay in the shape;
  // smaller alpha follows the point cloud's concavities more tightly, larger alpha tends to the convex hull.
  AlphaShape(double alpha) : alpha_(alpha) {
  }

  void addRoutes(const std::vector<MatchedRoute>& matched_routes) {
    constexpr Coordinate START_COORD{55.59784, 11.97298};
    constexpr double TOLERANCE_KM = 0.2;
    for (const auto& matched_route : matched_routes) {
      /*
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
      */
      matched_routes_.push_back(matched_route);

      auto simplified_route = simplifyCoordinates(matched_route.route.route);

      for (const auto& coord : simplified_route) {
        if (!coord.isValid()) {
          continue;
        }
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
    float min_x = std::numeric_limits<float>::max();
    float min_y = std::numeric_limits<float>::max();
    float max_x = std::numeric_limits<float>::lowest();
    float max_y = std::numeric_limits<float>::lowest();
    for (const auto& [lat, lon] : points_) {
      auto [x, y] = latLonToMeters(lat, lon);
      float fx = static_cast<float>(x);
      float fy = static_cast<float>(y);
      min_x = std::min(min_x, fx);
      min_y = std::min(min_y, fy);
      max_x = std::max(max_x, fx);
      max_y = std::max(max_y, fy);
      projected.emplace_back(fx, fy);
    }

    // Subdiv2D requires a rect that strictly contains every inserted point.
    cv::Rect2f bounds(static_cast<float>(min_x - 1.0), static_cast<float>(min_y - 1.0),
                      static_cast<float>(max_x - min_x + 2.0), static_cast<float>(max_y - min_y + 2.0));
    cv::Subdiv2D subdiv(bounds);
    for (const auto& p : projected) {
      try {
        subdiv.insert(p);
      } catch (const cv::Exception& e) {
        std::cerr << "Error inserting point into Subdiv2D: " << e.what() << std::endl;
        std::cerr << "Point causing error: (" << p.x << ", " << p.y << ")" << std::endl;
      }
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

  void paintMatches(const std::vector<WaySegment>& way_segments) {
    int count = 0;
    for (const auto& segment : way_segments) {
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
        double age = segment.getLastVisitAge();  // days
        // double thickness = std::clamp(5.0 * (1.0 - age / 365.0), 1.0, 5.0);  // Clamp thickness between 1 and 5
        // pixels
        double thickness = 2.0;  // Default thickness for the line
        cv::line(image_data_, cv::Point(x1, y1), cv::Point(x2, y2), segment.getLastVisitColor(),
                 static_cast<int>(thickness), cv::LINE_AA);
      }
    }
    // fprintf(stderr, "Painted %d way segments on tile z=%d, x=%d, y=%d\n", count, z_, x_, y_);
  }

  void paintHeatmap(const std::vector<MatchedRoute>& matched_routes) {
    // Traditional Strava heatmap gradient: dark red for lightly traveled pixels, through
    // orange and yellow, up to a white-hot core for the most heavily traveled ones.
    static const std::vector<cv::Scalar> heatmap_colors = {
        cv::Scalar(0, 0, 139, 255),  // dark red
        // cv::Scalar(0, 0, 255, 255),     // red
        cv::Scalar(0, 140, 255, 255),   // orange
        cv::Scalar(0, 255, 255, 255),   // yellow
        cv::Scalar(255, 255, 255, 255)  // white hot core
    };

    // Accumulate per-route coverage in a float buffer so pixels crossed by many different
    // activities build up brightness, then colorize with the orange-to-white-hot gradient.
    cv::Mat accumulator(256, 256, CV_32FC1, cv::Scalar(0));
    cv::Mat route_mask(256, 256, CV_8UC1);

    for (const auto& matched_route : matched_routes) {
      const auto& route = matched_route.route.route;
      if (route.size() < 2) {
        continue;
      }
      route_mask.setTo(cv::Scalar(0));
      bool drew_any = false;
      for (size_t i = 1; i < route.size(); ++i) {
        const auto& p1 = route[i - 1];
        const auto& p2 = route[i];
        int x1 = static_cast<int>((p1.lon - min_lon_) / (max_lon_ - min_lon_) * 256);
        int y1 = static_cast<int>((max_lat_ - p1.lat) / (max_lat_ - min_lat_) * 256);
        int x2 = static_cast<int>((p2.lon - min_lon_) / (max_lon_ - min_lon_) * 256);
        int y2 = static_cast<int>((max_lat_ - p2.lat) / (max_lat_ - min_lat_) * 256);
        // Skip segments that clearly miss this tile; cv::line clips the rest for us.
        if ((x1 < -8 && x2 < -8) || (x1 > 264 && x2 > 264) || (y1 < -8 && y2 < -8) || (y1 > 264 && y2 > 264)) {
          continue;
        }
        cv::line(route_mask, cv::Point(x1, y1), cv::Point(x2, y2), cv::Scalar(255), 1, cv::LINE_AA);
        drew_any = true;
      }
      if (drew_any) {
        cv::Mat route_mask_f;
        route_mask.convertTo(route_mask_f, CV_32FC1, 1.0 / 255.0);
        accumulator += route_mask_f;  // one contribution per route, so overlaps between activities stack
      }
    }

    // A small blur gives tracks the soft glow Strava's heatmap tiles have.
    int kernel_size = 1;
    if (z_ > 15) {
      kernel_size = 3;
    }
    if (kernel_size > 1) {
      cv::GaussianBlur(accumulator, accumulator, cv::Size(kernel_size, kernel_size), 0);
    }

    constexpr double max_value = 50;
    // cv::minMaxLoc(accumulator, nullptr, &max_value);
    // if (max_value <= 0.0) {
    //   return;
    // }

    for (int y = 0; y < accumulator.rows; ++y) {
      for (int x = 0; x < accumulator.cols; ++x) {
        float value = accumulator.at<float>(y, x);
        if (value <= 0.0f) {
          continue;
        }
        // Log scale keeps a single pass visible while heavily-traveled pixels saturate towards white.
        double normalized = std::clamp(std::log1p(value) / std::log1p(max_value), 0.0, 1.0);
        cv::Scalar color = color_map(normalized, heatmap_colors);
        auto& pixel = image_data_.at<cv::Vec4b>(y, x);
        pixel[0] = static_cast<uchar>(color[0]);
        pixel[1] = static_cast<uchar>(color[1]);
        pixel[2] = static_cast<uchar>(color[2]);
        pixel[3] = static_cast<uchar>(255);
      }
    }
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

  // Draw the full squadrat hex grid, regardless of which tiles were visited
  void paintGrid(double tile_size) {
    if (z_ < 11) {
      return;  // Only draw grid for zoom levels 11 and above
    }
    int alpha = 100;  // Adjust the alpha value for desired opacity
    if (z_ >= 14) {
      alpha = 200;  // Make grid lines more visible at higher zoom levels
    }
    cv::Scalar grey(20, 20, 20, alpha);
    for (const auto& [q, r] : hexesOverlappingTile(tile_size)) {
      auto polygon = hexToPixelPolygon(SquadratTile::hexVertices(q, r, tile_size));
      cv::polylines(image_data_, polygon, /*isClosed=*/true, grey, 1, cv::LINE_AA);
    }
  }

  void paint(const std::set<SquadratTile>& squadrat_tiles) {
    for (const auto& tile : squadrat_tiles) {
      auto polygon = hexToPixelPolygon(tile.getVertices());
      cv::fillConvexPoly(image_data_, polygon, tile.getLastVisitColor(), cv::LINE_AA);
      // Print the age at the hex's centroid
      cv::Point center(0, 0);
      for (const auto& p : polygon) {
        center += p;
      }
      center.x = center.x / static_cast<int>(polygon.size()) - 10;  // Center + shift text to be centered
      center.y = center.y / static_cast<int>(polygon.size()) + 5;
      if (z_ >= 11) {
        cv::putText(image_data_, std::to_string((int)tile.getLastVisitAge()) + "d", center, cv::FONT_HERSHEY_SIMPLEX,
                    0.4, cv::Scalar(0, 0, 0, 255), 1, cv::LINE_AA);
      }
    }
  }

  void paint(const AlphaShape& alpha_shape) {
    // paintRoute(alpha_shape.getMatchedRoutes());
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
  // Project hexagon corners (lat/lon) to pixel coordinates within this tile's 256x256 image
  std::vector<cv::Point> hexToPixelPolygon(const std::vector<Coordinate>& vertices) const {
    std::vector<cv::Point> polygon;
    polygon.reserve(vertices.size());
    for (const auto& v : vertices) {
      int px = static_cast<int>((v.lon - min_lon_) / (max_lon_ - min_lon_) * 256);
      int py = static_cast<int>((max_lat_ - v.lat) / (max_lat_ - min_lat_) * 256);
      polygon.emplace_back(px, py);
    }
    return polygon;
  }

  // Axial (q, r) coordinates of every hexagon whose bounding area could overlap this tile
  std::vector<std::pair<int, int>> hexesOverlappingTile(double tile_size) const {
    auto [min_x, min_y] = latLonToMeters(min_lat_, min_lon_);
    auto [max_x, max_y] = latLonToMeters(max_lat_, max_lon_);
    auto toAxial = [&](double x, double y) {
      double q = (HEX_SQRT3 / 3.0 * x - y / 3.0) / tile_size;
      double r = (2.0 / 3.0 * y) / tile_size;
      return std::make_pair(q, r);
    };
    auto [q1, r1] = toAxial(min_x, min_y);
    auto [q2, r2] = toAxial(min_x, max_y);
    auto [q3, r3] = toAxial(max_x, min_y);
    auto [q4, r4] = toAxial(max_x, max_y);
    // Pad by one extra ring so hexagons whose center falls outside the tile still get drawn.
    int q_start = static_cast<int>(std::floor(std::min({q1, q2, q3, q4}))) - 1;
    int q_end = static_cast<int>(std::ceil(std::max({q1, q2, q3, q4}))) + 1;
    int r_start = static_cast<int>(std::floor(std::min({r1, r2, r3, r4}))) - 1;
    int r_end = static_cast<int>(std::ceil(std::max({r1, r2, r3, r4}))) + 1;

    std::vector<std::pair<int, int>> hexes;
    for (int q = q_start; q <= q_end; ++q) {
      for (int r = r_start; r <= r_end; ++r) {
        hexes.emplace_back(q, r);
      }
    }
    return hexes;
  }

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
    route.activity_id = feature["properties"]["activity_id"];
    if (activity_id_blacklist.end() !=
        std::find(activity_id_blacklist.begin(), activity_id_blacklist.end(), route.activity_id)) {
      continue;
    }
    route.name = feature["properties"]["name"];
    route.date = feature["properties"]["date"];
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

// Decode polyline geometry
static std::vector<Coordinate> decodePolyline(const std::string& encoded, double precision) {
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
    coordinates.push_back(Coordinate{lat / precision, lng / precision});
  }

  return coordinates;
}

class RouteMatcher {
 public:
  RouteMatcher() {
    char* valhalla_config_file = std::getenv("VALHALLA_CONFIG_FILE");
    std::string valhalla_config_file_str(valhalla_config_file ? valhalla_config_file : VALHALLA_CONFIG_FILE);
    fprintf(stderr, "Initializing Valhalla with config from: %s\n", valhalla_config_file_str.c_str());
    TimerLog timer("Valhalla initialization");
    const auto& config = valhalla::config(valhalla_config_file_str);
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
    double visit_time = route.getMilliseconds();
    for (const auto& trip_route : api.trip().routes()) {
      for (const auto& leg : trip_route.legs()) {
        auto leg_shape = decodePolyline(leg.shape(), 1e6);
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
          matched_route.way_segments.back().visit(visit_time);
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

  void clearCache() {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    tile_cache_.clear();
  }

  virtual void addRoutes(const std::vector<MatchedRoute>& matched_routes) = 0;

 private:
  std::mutex cache_mutex_;
  std::unordered_map<TileKey, Tile, TileKey::Hash> tile_cache_;
};

class AlphaShapeTileGenerator : public TileGenerator {
 public:
  AlphaShapeTileGenerator(const std::vector<MatchedRoute>& matched_routes, double alpha) : alpha_shape_(alpha) {
    addRoutes(matched_routes);
  }

  void addRoutes(const std::vector<MatchedRoute>& matched_routes) override {
    std::scoped_lock<std::shared_mutex> full_lock(shared_mutex_);
    alpha_shape_.addRoutes(matched_routes);
  }

  Tile generateTile(int z, int x, int y) override {
    Tile tile(z, x, y);
    {
      std::shared_lock<std::shared_mutex> shared_lock(shared_mutex_);
      tile.paint(alpha_shape_);
    }
    return tile;
  }

 private:
  std::shared_mutex shared_mutex_;
  AlphaShape alpha_shape_;
};

class SquadratTileGenerator : public TileGenerator {
 public:
  SquadratTileGenerator(const std::vector<MatchedRoute>& matched_routes, double tile_size_raw)
      : tile_size_(meter2size(tile_size_raw)) {
    addRoutes(matched_routes);
  }

  void addRoutes(const std::vector<MatchedRoute>& matched_routes) override {
    for (const auto& matched_route : matched_routes) {
      double visit_time = matched_route.route.getMilliseconds();
      for (const auto& coordinate : matched_route.route.route) {
        std::scoped_lock<std::shared_mutex> full_lock(shared_mutex_);
        auto it = squadrat_tiles_.emplace(coordinate.lat, coordinate.lon, tile_size_);
        it.first->visit(visit_time);
      }
    }
    clearCache();
  }

  Tile generateTile(int z, int x, int y) override {
    Tile tile(z, x, y);
    {
      std::shared_lock<std::shared_mutex> shared_lock(shared_mutex_);
      tile.paint(squadrat_tiles_);
    }
    tile.paintGrid(tile_size_);
    return tile;
  }

 private:
  std::shared_mutex shared_mutex_;
  std::set<SquadratTile> squadrat_tiles_;
  double tile_size_;
};

class TraversalTileGenerator : public TileGenerator {
 public:
  TraversalTileGenerator(const std::vector<MatchedRoute>& matched_routes) {
    addRoutes(matched_routes);
  }

  void addRoutes(const std::vector<MatchedRoute>& matched_routes) override {
    std::unordered_map<std::size_t, WaySegment> traversal_counts;
    for (const auto& matched_route : matched_routes) {
      double visit_time = matched_route.route.getMilliseconds();
      for (const auto& segment : matched_route.way_segments) {
        auto it = traversal_counts.find(segment.edge_id);
        if (it == traversal_counts.end()) {
          traversal_counts[segment.edge_id] = segment;
          traversal_counts[segment.edge_id].visit(visit_time);
        } else {
          it->second.traversal_count += 1;
          it->second.visit(visit_time);
        }
      }
    }
    // Sort after visit time and insert in traversal_counts_
    std::scoped_lock<std::shared_mutex> full_lock(shared_mutex_);
    for (const auto& [key, segment] : traversal_counts) {
      traversal_counts_.push_back(segment);
    }
    std::sort(traversal_counts_.begin(), traversal_counts_.end(),
              [](const WaySegment& a, const WaySegment& b) { return a.getLastVisitTime() < b.getLastVisitTime(); });
    clearCache();
  }

  Tile generateTile(int z, int x, int y) override {
    Tile tile(z, x, y);
    {
      std::shared_lock<std::shared_mutex> shared_lock(shared_mutex_);
      tile.paintMatches(traversal_counts_);
    }
    return tile;
  }

 private:
  std::shared_mutex shared_mutex_;
  std::vector<WaySegment> traversal_counts_;
};

class StravaHeatmapTileGenerator : public TileGenerator {
 public:
  StravaHeatmapTileGenerator(const std::vector<MatchedRoute>& matched_routes) {
    addRoutes(matched_routes);
  }

  void addRoutes(const std::vector<MatchedRoute>& matched_routes) override {
    std::scoped_lock<std::shared_mutex> full_lock(shared_mutex_);
    matched_routes_.insert(matched_routes_.end(), matched_routes.begin(), matched_routes.end());
    clearCache();
  }

  Tile generateTile(int z, int x, int y) override {
    Tile tile(z, x, y);
    {
      std::shared_lock<std::shared_mutex> shared_lock(shared_mutex_);
      tile.paintHeatmap(matched_routes_);
    }
    return tile;
  }

 private:
  std::shared_mutex shared_mutex_;
  std::vector<MatchedRoute> matched_routes_;
};

static void persistHeatmap(const std::vector<MatchedRoute>& matched_routes, const std::string& heatmap_file_path_str) {
  witxheatmap::PHeatmap heatmap_proto;
  for (const auto& matched_route : matched_routes) {
    *heatmap_proto.add_matched_routes() = matchedRouteToProto(matched_route);
  }
  std::ofstream heatmap_file(heatmap_file_path_str, std::ios::binary);
  if (!heatmap_file.is_open()) {
    throw std::runtime_error("Failed to open heatmap file for writing: " + heatmap_file_path_str);
  }
  if (!heatmap_proto.SerializeToOstream(&heatmap_file)) {
    throw std::runtime_error("Failed to serialize heatmap to: " + heatmap_file_path_str);
  }
  heatmap_file.close();
}

class ActivityIngester {
 public:
  ActivityIngester(const std::vector<TileGenerator*>& tile_generators) : tile_generators_(tile_generators) {
    char* ingest_folder = std::getenv("INGEST_FOLDER");
    ingest_folder_str_ = ingest_folder ? ingest_folder : DEFAULT_INGEST_FOLDER;
    // On creation re-ingest the already ingested activities
    ingestFolder(ingest_folder_str_ + "/ingested", false);
    start();
  }

  ~ActivityIngester() {
    if (ingest_thread_.joinable()) {
      ingest_thread_.join();
    }
  }

  void ingestFolder(const std::string& ingest_folder_str, bool move_activities) {
    if (!std::filesystem::exists(ingest_folder_str)) {
      std::cerr << "Ingest folder does not exist: " << ingest_folder_str << std::endl;
      return;
    }
    std::filesystem::directory_iterator ingest_dir(ingest_folder_str);

    std::vector<Route> routes;
    for (const auto& entry : ingest_dir) {
      // Found new json file
      if (entry.is_regular_file() && entry.path().extension() == ".json") {
        std::cout << "Found new file: " << entry.path() << std::endl;
        nlohmann::json activity_json = nlohmann::json::parse(std::ifstream(entry.path()));
        Route route;
        route.name = activity_json.value("name", "");
        route.date = activity_json.value("start_date", "");
        route.date_format = 1;
        std::string polyline = activity_json.at("map").value("polyline", "");
        route.route = decodePolyline(polyline, 1e5);
        std::cout << "Ingesting activity: " << route.name << " at time " << route.date << " with " << route.route.size()
                  << " points" << std::endl;
        routes.push_back(route);

        // Move to ingested folder
        if (move_activities) {
          std::filesystem::path ingested_folder = ingest_folder_str + "/ingested";
          if (!std::filesystem::exists(ingested_folder)) {
            std::filesystem::create_directory(ingested_folder);
          }
          std::filesystem::rename(entry.path(), ingested_folder / entry.path().filename());
        }
      }
    }
    if (!routes.empty()) {
      std::vector<MatchedRoute> matches_routes = ingest_matcher_.matchAllRoutes(routes);
      for (const auto& tile_generator : tile_generators_) {
        tile_generator->addRoutes(matches_routes);
      }
    }
  }

  void start() {
    ingest_thread_ = std::thread([this]() {
      std::cout << "Ingest folder: " << ingest_folder_str_ << std::endl;

      while (true) {
        // Clear cache once a day
        if (millisecondsNow() - last_cleared_time_ > MILLISECONDS_PER_DAY) {
          // Clear cache logic here
          for (const auto& tile_generator : tile_generators_) {
            tile_generator->clearCache();
          }
          last_cleared_time_ = millisecondsNow();
        }

        ingestFolder(ingest_folder_str_, true);
        std::this_thread::sleep_for(std::chrono::seconds(60));
      }
    });
  }

 private:
  RouteMatcher ingest_matcher_;
  std::string ingest_folder_str_;
  std::thread ingest_thread_;
  double last_cleared_time_ = millisecondsNow();
  std::vector<TileGenerator*> tile_generators_;
};

int main(int argc, char** argv) {
  std::string url_path = "/tiles";
  int port = 9090;
  if (argc > 2) {
    url_path = argv[1];
    port = std::stoi(argv[2]);
  }

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

  char* heatmap_file_path = std::getenv("HEATMAP_FILE_PATH");
  std::string heatmap_file_path_str(heatmap_file_path ? heatmap_file_path : HEATMAP_FILE);

  std::vector<MatchedRoute> matched_routes;
  {
#if !RESTORE
    auto routes = load_all_routes();
    RouteMatcher matcher;
    matched_routes = matcher.matchAllRoutes(routes);
    persistHeatmap(matched_routes, heatmap_file_path_str);

#else

    std::ifstream heatmap_file(heatmap_file_path_str, std::ios::binary);
    if (!heatmap_file.is_open()) {
      throw std::runtime_error("Failed to open heatmap file for reading: " + heatmap_file_path_str);
    }
    TimerLog restore_timer("Loading matched routes from heatmap file");
    witxheatmap::PHeatmap heatmap_proto;
    if (!heatmap_proto.ParseFromIstream(&heatmap_file)) {
      throw std::runtime_error("Failed to parse heatmap file: " + heatmap_file_path_str);
    }
    matched_routes.reserve(heatmap_proto.matched_routes_size());
    for (const auto& proto_matched_route : heatmap_proto.matched_routes()) {
      matched_routes.push_back(matchedRouteFromProto(proto_matched_route));
    }
#endif
  }
  std::vector<TileGenerator*> tile_generators;

  std::unordered_map<int, std::unique_ptr<AlphaShapeTileGenerator>> alpha_shapes;
  {
    TimerLog alpha_shapes_timer("Generating alpha shapes for radius 4000");
    auto it = alpha_shapes.emplace(
        4000, std::make_unique<AlphaShapeTileGenerator>(matched_routes, static_cast<double>(4000)));
    tile_generators.push_back(it.first->second.get());
  }
  {
    TimerLog alpha_shapes_timer("Generating alpha shapes for radius 7000");
    auto it = alpha_shapes.emplace(
        7000, std::make_unique<AlphaShapeTileGenerator>(matched_routes, static_cast<double>(7000)));
    tile_generators.push_back(it.first->second.get());
  }
  {
    TimerLog alpha_shapes_timer("Generating alpha shapes for radius 10000");
    auto it = alpha_shapes.emplace(
        10000, std::make_unique<AlphaShapeTileGenerator>(matched_routes, static_cast<double>(10000)));
    tile_generators.push_back(it.first->second.get());
  }
  std::mutex alpha_shapes_mutex;

  // Main route planning endpoint
  svr.Get(url_path + R"(/coverage/(\d+)/(\d+)/(\d+).png)", [&matched_routes, &alpha_shapes, &alpha_shapes_mutex](
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
    TimerLog squadrat_tiles_timer("Generating squadrats for radius 500");
    auto it = squadrat_tile_generators.emplace(
        500, std::make_unique<SquadratTileGenerator>(matched_routes, static_cast<double>(500)));
    tile_generators.push_back(it.first->second.get());
  }
  {
    TimerLog squadrat_tiles_timer("Generating squadrats for radius 1000");
    auto it = squadrat_tile_generators.emplace(
        1000, std::make_unique<SquadratTileGenerator>(matched_routes, static_cast<double>(1000)));
    tile_generators.push_back(it.first->second.get());
  }
  {
    TimerLog squadrat_tiles_timer("Generating squadrats for radius 1600");
    auto it = squadrat_tile_generators.emplace(
        1600, std::make_unique<SquadratTileGenerator>(matched_routes, static_cast<double>(1600)));
    tile_generators.push_back(it.first->second.get());
  }
  {
    TimerLog squadrat_tiles_timer("Generating squadrats for radius 2000");
    auto it = squadrat_tile_generators.emplace(
        2000, std::make_unique<SquadratTileGenerator>(matched_routes, static_cast<double>(2000)));
    tile_generators.push_back(it.first->second.get());
  }
  std::mutex squadrat_tiles_mutex;

  svr.Get(url_path + R"(/squadrat/(\d+)/(\d+)/(\d+).png)", [&matched_routes, &squadrat_tile_generators,
                                                            &squadrat_tiles_mutex](const httplib::Request& req,
                                                                                   httplib::Response& res) {
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

  TraversalTileGenerator traversal_tile_generator(matched_routes);
  tile_generators.push_back(&traversal_tile_generator);

  svr.Get(url_path + R"(/traversal/(\d+)/(\d+)/(\d+).png)",
          [&traversal_tile_generator](const httplib::Request& req, httplib::Response& res) {
            // Heatmap XYZ tile request from url like /traversal/{z}/{x}/{y}.png
            for (const auto& param : req.path_params) {
              fprintf(stderr, "Path param: %s = %s\n", param.first.c_str(), param.second.c_str());
            }
            int z = std::stoi(req.matches[1]);
            int x = std::stoi(req.matches[2]);
            int y = std::stoi(req.matches[3]);
            std::string user = req.has_param("user") ? req.get_param_value("user") : "";

            // For now just return a placeholder PNG image (256x256 pixel)
            Tile tile = traversal_tile_generator.getTile(z, x, y);
            cv::Mat png_data = tile.getImage();
            std::vector<unsigned char> buffer;
            cv::imencode(".png", png_data, buffer);
            res.set_content(reinterpret_cast<const char*>(buffer.data()), buffer.size(), "image/png");
          });

  StravaHeatmapTileGenerator heatmap_tile_generator(matched_routes);
  tile_generators.push_back(&heatmap_tile_generator);

  svr.Get(url_path + R"(/heatmap/(\d+)/(\d+)/(\d+).png)",
          [&heatmap_tile_generator](const httplib::Request& req, httplib::Response& res) {
            // Heatmap XYZ tile request from url like /heatmap/{z}/{x}/{y}.png
            for (const auto& param : req.path_params) {
              fprintf(stderr, "Path param: %s = %s\n", param.first.c_str(), param.second.c_str());
            }
            int z = std::stoi(req.matches[1]);
            int x = std::stoi(req.matches[2]);
            int y = std::stoi(req.matches[3]);
            std::string user = req.has_param("user") ? req.get_param_value("user") : "";

            // For now just return a placeholder PNG image (256x256 pixel)
            Tile tile = heatmap_tile_generator.getTile(z, x, y);
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

  ActivityIngester ingester(tile_generators);

  // Start server
  std::cout << "\n==================================================" << std::endl;
  std::cout << "Route Planning API Server" << std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << "Server starting on http://0.0.0.0:" << port << std::endl;
  std::cout << "  GET /health - Health check" << std::endl;
  std::cout << "  GET /tile?z=0&x=0&y=0&radius=10" << std::endl;
  std::cout << "\nParameters:" << std::endl;
  std::cout << "  z       - Zoom level (default: 0)" << std::endl;
  std::cout << "  x       - Tile x coordinate (default: 0)" << std::endl;
  std::cout << "  y       - Tile y coordinate (default: 0)" << std::endl;
  std::cout << "  radius  - Search radius in km (default: 10)" << std::endl;

  svr.listen("0.0.0.0", port);

  return 0;
}