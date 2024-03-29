cmake_minimum_required(VERSION 3.0)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 11)

if(WIN32)
	set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${3RDPARTY_PATH}/zlib ${3RDPARTY_PATH}/libjpeg ${3RDPARTY_PATH}/libpng ${3RDPARTY_PATH}/openexr)
	set(JPEG_NAMES "libjpeg")
endif()

find_package(PNG)
find_package(JPEG)
find_package(TIFF)
find_package(OpenEXR)

add_library(image STATIC include/image.h src/image_io.cc src/image_manip.cc)
target_compile_options(image PRIVATE "-w")
target_include_directories(image PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
if(PNG_FOUND)
	target_compile_definitions(image PUBLIC -DWITH_PNG)
	target_include_directories(image PRIVATE ${PNG_INCLUDE_DIRS})
	target_link_libraries(image ${PNG_LIBRARIES})
endif()
if(JPEG_FOUND)
	target_compile_definitions(image PUBLIC -DWITH_JPEG)
	target_include_directories(image PRIVATE ${JPEG_INCLUDE_DIR})
	target_link_libraries(image ${JPEG_LIBRARIES})
endif()
if(TIFF_FOUND)
	target_compile_definitions(image PUBLIC -DWITH_TIFF)
	target_include_directories(image PRIVATE ${TIFF_INCLUDE_DIR})
	target_link_libraries(image ${TIFF_LIBRARIES})
endif()
#if(OPENEXR_FOUND)
#	target_compile_definitions(image PUBLIC -DWITH_EXR)
#	target_include_directories(image PRIVATE ${OPENEXR_INCLUDE_DIRS})
#	target_link_libraries(image ${OPENEXR_LIBRARIES})
#endif()
