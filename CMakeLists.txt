cmake_minimum_required(VERSION 3.10)
project(gravel_heatmap)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-stringop-overflow")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wno-stringop-overflow")

# Find OpenMP
find_package(OpenMP REQUIRED)

# Find OSRM package
find_package(PkgConfig REQUIRED)
pkg_check_modules(OSRM REQUIRED libosrm)

# Find Boost (required by OSRM)
find_package(Boost REQUIRED COMPONENTS filesystem system iostreams thread)

# Add executable
add_executable(witx-heatmap witx-heatmap.cpp)

# Include directories
target_include_directories(witx-heatmap PRIVATE 
    ${OSRM_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(witx-heatmap 
    ${OSRM_LIBRARIES}
    ${Boost_LIBRARIES}
    pthread
    OpenMP::OpenMP_CXX
)

# Compiler flags
target_compile_options(witx-heatmap PRIVATE ${OSRM_CFLAGS_OTHER})

