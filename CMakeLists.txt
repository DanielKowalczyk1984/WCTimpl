CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(WCTproject)

MESSAGE("${PMCproject_BINARY_DIR}")

set(FindFiles ${WCTproject_DIR}/FindFiles)

include(FindFiles/FindGlib.cmake)

include(FindFiles/FindGurobi.cmake)

find_package(Boost 1.60.0)

include_directories(${GUROBI_INCLUDE_DIR} ${GLIB_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

SET(CMAKE_C_FLAGS "-g -std=c99  -D_XOPEN_SOURCE=500 -pedantic -Wall -Wshadow -W -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wpointer-arith -Wnested-externs -Wundef -Wcast-qual -Wcast-align -Wwrite-strings")

SET(CMAKE_CXX_FLAGS "-g -std=c++11")

FILE(GLOB SOURCE  "src/*.c" "src/*.cc")

set(GLOB headers "src/*.h" ${GUROBI_INCLUDE_DIR} "src/*.hpp")

add_executable(wct ${SOURCE} ${headers})

target_link_libraries(wct ${GUROBI_LIBRARY} ${GLIB_LIBRARIES} ${Boost_LIBRARIES} -lm -lpthread -lgurobi60)
