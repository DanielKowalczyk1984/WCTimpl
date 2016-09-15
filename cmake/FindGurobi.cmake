if (GUROBI_INCLUDE_DIR)
  # in cache already
  set(GUROBI_FOUND TRUE)
  set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
  set(GUROBI_LIBRARIES "${GUROBI_LIBRARY};${GUROBI_CXX_LIBRARY}" )
else (GUROBI_INCLUDE_DIR)

find_path(GUROBI_INCLUDE_DIR
          NAMES  gurobi_c++.h
          PATHS   
          "/opt/gurobi652/linux64/include"
          "/opt/gurobi650/linux64/include/"
                "$ENV{GUROBI_HOME}/include"
                  "/Library/gurobi502/mac64/include"
                  "/Library/gurobi650/mac64/include/"
                  "/Library/gurobi604/mac64/include/"
                 "C:\\libs\\gurobi502\\include"
                 "/opt/gurobi563/linux64/include"
          )

find_library( GUROBI_LIBRARY
              NAMES
              gurobi65
              gurobi60
              gurobi
              gurobi45
              gurobi46
              gurobi50
              gurobi51
              gurobi52
              gurobi55
              gurobi56
              gurobi60
              gurobi563
              PATHS
              "/opt/gurobi652/linux64/lib/"
              "/opt/gurobi650/linux64/lib/"
                    "/Library/gurobi604/mac64/lib/"
                    "/opt/gurobi600/linux64/lib/"
                    "$ENV{GUROBI_HOME}/lib"
                    "/Library/gurobi650/mac64/lib"
                    "/Library/gurobi502/mac64/lib"
                    "C:\\libs\\gurobi502\\lib"
                    "/opt/gurobi563/linux64/lib"
              )

find_library( GUROBI_CXX_LIBRARY
              NAMES
              gurobi_c++
              libgurobi56
              libgurobi

              PATHS
              "/opt/gurobi652/linux64/lib/"
              "/opt/gurobi650/linux64/lib/"
                    "$ENV{GUROBI_HOME}/lib"
                    "/Library/gurobi650/mac64/lib"
                    "/Library/gurobi604/mac64/lib"
                    "/Library/gurobi502/mac64/lib"
                    "C:\\libs\\gurobi502\\lib"
                    "/opt/gurobi563/linux64/lib/"
              )

set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY}" )

message(${GUROBI_CXX_LIBRARY})
message(${GUROBI_INCLUDE_DIR})
message(${GUROBI_LIBRARY})
# use c++ headers as default
# set(GUROBI_COMPILER_FLAGS "-DIL_STD" CACHE STRING "Gurobi Compiler Flags")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
# if all listed variables are TRUE
# find_package_handle_standard_args(GUROBI DEFAULT_MSG
#                                    GUROBI_CXX_LIBRARY GUROBI_INCLUDE_DIR)

mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY GUROBI_CXX_LIBRARY)
set (CMAKE_SHARED_LINKER_FLAGS "-lgurobi_c++ -lgurobi65, --as-needed")

endif(GUROBI_INCLUDE_DIR)