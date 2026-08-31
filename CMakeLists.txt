cmake_minimum_required(VERSION 3.21)
project(gravel_heatmap)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-stringop-overflow")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wno-stringop-overflow")

option(ENABLE_ASAN "Build with AddressSanitizer" OFF)
if(ENABLE_ASAN)
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address)
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)
find_package(OpenCV REQUIRED)

# Valhalla's ptree config API only needs Boost headers (property_tree).
find_package(Boost 1.71 REQUIRED)

# Use the prebuilt Valhalla install (submodule) for map matching instead of OSRM.
set(VALHALLA_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/valhalla/install)
set(ENV{PKG_CONFIG_PATH} "${VALHALLA_INSTALL_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")

find_package(PkgConfig REQUIRED)
pkg_check_modules(VALHALLA REQUIRED libvalhalla)
# libvalhalla is a static archive; its Requires.private aren't pulled in by
# pkg_check_modules unless using --static, which over-recurses into every
# dependency's own private/static-only deps (e.g. curl's bundled ssh/ldap/etc
# libs that have no dev symlinks installed here). Query the direct private
# deps by name instead, using normal (dynamic) linking for all of them.
pkg_check_modules(VALHALLA_DEPS REQUIRED
    spatialite
    sqlite3
    luajit
    geos
    openssl
    libcurl
    libtiff-4
    liblz4
)


find_path(PROTOZERO_INCLUDE_DIR protozero/pbf_reader.hpp)

# Add executable
add_executable(witx-heatmap witx-heatmap.cpp)

target_include_directories(witx-heatmap PRIVATE
    ${VALHALLA_INCLUDE_DIRS}
    ${VALHALLA_DEPS_INCLUDE_DIRS}
    ${VALHALLA_STATIC_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
)
target_compile_options(witx-heatmap PRIVATE ${VALHALLA_CFLAGS_OTHER} ${VALHALLA_DEPS_CFLAGS_OTHER} ${VALHALLA_STATIC_CFLAGS_OTHER})

# Link libraries
target_link_directories(witx-heatmap PRIVATE ${VALHALLA_LIBRARY_DIRS} ${VALHALLA_DEPS_LIBRARY_DIRS})
target_link_libraries(witx-heatmap
    ${VALHALLA_LIBRARIES}
    ${VALHALLA_DEPS_LIBRARIES}
    geotiff
    Boost::boost
    pthread
    OpenMP::OpenMP_CXX
    ${OpenCV_LIBS}
)

