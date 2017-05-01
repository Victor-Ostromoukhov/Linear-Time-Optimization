cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 OLD)

##### PROJECT #####
project(discrepancy)

##### OPTIONS #####
set(EXECUTABLE_OUTPUT_PATH ./bin)
find_package(Boost 1.42.0 COMPONENTS program_options REQUIRED)
#find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS
    ${CMAKE_CXX_FLAGS}
    # ${OpenMP_CXX_FLAGS}
)

include_directories(
	./src
	${Boost_INCLUDE_DIR}
)

link_directories(
	.
	${Boost_LIBRARY_DIRS}
)

add_executable(discrepancy src/parallelism.cpp src/Points_2D.cpp src/Points_1D.cpp src/Points_3D.cpp)
target_link_libraries(
	discrepancy
	${Boost_LIBRARIES}
)

