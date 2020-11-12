cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

project(ReconComp)
file(GLOB src
    "slic.cc"
    "raw.cc"
    "main.cc"
)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11 -lstdc++fs -w -funroll-loops -ffast-math -fassociative-math -funsafe-math-optimizations -ftree-vectorize") #-march=native -O2 -flto

add_subdirectory(image)
find_package(OpenCV REQUIRED)
find_package(LibRaw REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP)
find_package(Threads)
if(OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
add_executable(rc ${src})
target_link_libraries(rc image ${OpenCV_LIBS} ${LibRaw_LIBRARIES} Eigen3::Eigen)
target_compile_definitions(rc PRIVATE ${LibRaw_DEFINITIONS})
target_include_directories(rc PUBLIC ${LibRaw_INCLUDE_DIR} "/usr/include/opencv2" ${OpenCV_INCLUDE_DIRS} "${OpenCV_INCLUDE_DIRS}/opencv2")

