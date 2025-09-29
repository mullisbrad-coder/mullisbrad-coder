cmake_minimum_required(VERSION 3.15)
project(C3DConvexPoly LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(c3dconvexpoly STATIC
    src/C3DConvexPoly.cpp
)

target_include_directories(c3dconvexpoly PUBLIC include)

add_executable(convex_example examples/usage.cpp)
target_link_libraries(convex_example PRIVATE c3dconvexpoly)

