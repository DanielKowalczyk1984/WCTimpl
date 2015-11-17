CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(WCTproject)

MESSAGE("${PMCproject_BINARY_DIR}")

include(${WCTproject_BINARY_DIR}/FindGurobi.cmake)

include(${WCTproject_BINARY_DIR}/FindGlib.cmake)

include_directories(${GUROBI_INCLUDE_DIR} ${GLIB_INCLUDE_DIRS})

SET(CMAKE_C_FLAGS "-g -std=c99  -D_XOPEN_SOURCE=500 -pedantic -Wall -Wshadow -W -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wpointer-arith -Wnested-externs -Wundef -Wcast-qual -Wcast-align -Wwrite-strings")

SET(CMAKE_CXX_FLAGS "-g -std=c++11")

FILE(GLOB SOURCE  "src/*.c" "src/*.cc")

set(GLOB headers "src/*.h" ${GUROBI_INCLUDE_DIR} "src/*.hpp")

add_executable(wct ${SOURCE} ${headers})

target_link_libraries(wct ${GUROBI_LIBRARY} ${GLIB_LIBRARIES} -lm -lpthread -lgurobi60)
