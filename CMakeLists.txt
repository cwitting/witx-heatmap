cmake_minimum_required(VERSION 3.21)
project(gravel_heatmap)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-stringop-overflow")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wno-stringop-overflow")

# Find OpenMP
find_package(OpenMP REQUIRED)

# Build OSRM from the in-tree submodule instead of relying on a system
# install / pkg-config. This exposes the osrm::osrm target (and friends).
add_subdirectory(osrm-backend)

# Add executable
add_executable(witx-heatmap witx-heatmap.cpp)

# Include directories (OSRM headers + the version.hpp it generates at configure time)
target_include_directories(witx-heatmap PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/osrm-backend/include
    ${CMAKE_CURRENT_BINARY_DIR}/osrm-backend/include
    ${CMAKE_CURRENT_SOURCE_DIR}/osrm-backend/third_party/fmt/include
)

# OSRM's json_renderer.hpp uses fmt in header-only mode; match that here so we
# don't pull in the system libfmt (which requires linking a compiled fmt lib).
target_compile_definitions(witx-heatmap PRIVATE FMT_HEADER_ONLY)

# Link libraries (osrm::osrm pulls in its Boost/TBB/LibArchive deps transitively)
target_link_libraries(witx-heatmap
    osrm
    pthread
    OpenMP::OpenMP_CXX
)

