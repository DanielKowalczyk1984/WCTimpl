CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(WCTproject)

MESSAGE("${PMCproject_BINARY_DIR}")

SET(FindFiles ${WCTproject_DIR}/FindFiles)

include(FindFiles/FindGlib.cmake)

include(FindFiles/FindGurobi.cmake)

find_package(Boost)

include_directories(${GUROBI_INCLUDE_DIR} ${GLIB_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

SET(CMAKE_C_FLAGS " -std=c99 -O5 -D_XOPEN_SOURCE=500 -pedantic -Wall -Wshadow -W -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wpointer-arith -Wnested-externs -Wundef -Wcast-qual -Wcast-align -Wwrite-strings")

SET(CMAKE_CXX_FLAGS " -std=c++11 -O5")

FILE(GLOB SOURCE  "src/*.c" "src/*.cc")

SET(GLOB headers "src/*.h" ${GUROBI_INCLUDE_DIR} "src/*.hpp")

add_executable(wct ${SOURCE} ${headers})

target_link_libraries(wct ${GUROBI_LIBRARY} ${GLIB_LIBRARIES} ${Boost_LIBRARIES} -lm -lpthread -lgurobi60)
