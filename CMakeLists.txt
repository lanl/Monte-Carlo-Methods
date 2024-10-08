cmake_minimum_required(VERSION 3.8)

project(main CXX)

set(CMAKE_CXX_STANDARD 23)

add_subdirectory(extern/matplotlib-cpp)

add_executable(main 
  carlo/model/ising.cpp
  carlo/image.cpp
  main.cpp 
)
target_include_directories(main PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/gif-h
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/stb
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/matplotlib-cpp
)
target_compile_options(main PRIVATE -Wno-deprecated)
target_link_libraries(main PRIVATE matplotlib_cpp)