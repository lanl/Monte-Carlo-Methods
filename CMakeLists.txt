cmake_minimum_required(VERSION 3.8)

project(main CXX)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(extern/matplotlib-cpp)

add_executable(main 
  carlo/model/ising.hpp
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
target_link_libraries(main PRIVATE matplotlib_cpp fftw3)

# Because my python binder complains
target_compile_options(main PRIVATE -Wno-deprecated-declarations)