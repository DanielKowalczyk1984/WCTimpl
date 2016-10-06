CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(WCTproject)

MESSAGE("${PMCproject_BINARY_DIR}")

SET(cmake ${WCTproject_DIR}/cmake)

include(cmake/FindGlib.cmake)

include(cmake/FindGurobi.cmake)

find_package(Boost)

include_directories(${GUROBI_INCLUDE_DIR} ${GLIB_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -DNDEBUG -O5 -Wall -D_XOPEN_SOURCE=500")

SET(CMAKE_C_FLAGS_RELEASE "-O5 -Wall -D_XOPEN_SOURCE=500 -pedantic -std=c99")

SET(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -O5 -g -Wall -D_XOPEN_SOURCE=500 ")

SET(CMAKE_C_FLAGS_DEBUG "-O5 -g -Wall -D_XOPEN_SOURCE=500 -pedantic -std=c99")

FILE(GLOB SOURCE  "src/*.c" "src/*.cc")

SET(GLOB headers "src/*.h" ${GUROBI_INCLUDE_DIR} "src/*.hpp")

add_executable(wct ${SOURCE} ${headers})

target_link_libraries(wct ${GUROBI_LIBRARY} ${GLIB_LIBRARIES} ${Boost_LIBRARIES} -lm -lpthread)
