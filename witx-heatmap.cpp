// #include <omp.h>

#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <osrm/coordinate.hpp>
#include <osrm/engine_config.hpp>
#include <osrm/json_container.hpp>
#include <osrm/nearest_parameters.hpp>
#include <osrm/osrm.hpp>
#include <osrm/route_parameters.hpp>
#include <osrm/status.hpp>
#include <stdexcept>
#include <vector>

#include "json.hpp"

struct Coordinate {
  double lat{};
  double lon{};
};

using Route = std::vector<Coordinate>;

#define ROUTE_FILE "/home/christian/git/witx-heatmap/data/christian_strava.geojson"

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

// ============================================================================
// Main program
// ============================================================================
int main(int argc, char* argv[]) {
  // Print program name
  fprintf(stderr, "Witx Heatmap Route Planner\n");
  auto routes = load_all_routes();
  int route_count = 0;
  for (const auto& route : routes) {
    fprintf(stderr, "Route with %zu points\n", route.size());
    route_count++;
    fprintf(stderr, "Route count: %d\n", route_count);
  }
  fprintf(stderr, "Total routes: %d\n", route_count);
}
