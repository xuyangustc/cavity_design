cmake_minimum_required (VERSION 2.6)
project (migrated)

# The version number.
set (Dataio_VERSION_MAJOR 0)
set (Dataio_VERSION_MINOR 0)

# configure a header file to pass some of the CMake settings
# to the source code
#configure_file (
#  "${PROJECT_SOURCE_DIR}/dataioconfig.h.in"
#  "${PROJECT_BINARY_DIR}/dataioconfig.h"
#  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include/" "/home/xuyang/lib/")
SET(ENV{LD_LIBRARY_PATH} "${PROJECT_BINARY_DIR}/lib:/usr/local/lib:/usr/local/lib64:/home/xuyang/lib:$ENV{LD_LIBRARY_PATH}")
SET(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.66.0 REQUIRED COMPONENTS system filesystem timer)
add_definitions(-D_BOOST_NO_CXX11_SCOPED_ENUMS -D_Boost_USE_STATIC_LIBS=ON)
#SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -std=c++11 -lmlpack -larmadillo -lblas -lhdf5_cpp -lhdf5 -lboost_filesystem -lboost_system" )
#SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -std=c++11 -lgslcblas -lboost_filesystem -lboost_system -lboost_timer -lgsl -lloop" )
#SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -std=c++11 -L/usr/local/lib/ -lpthread -lgslcblas -lboost_filesystem -lboost_system -lboost_timer -lgsl -lloop" )
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O4 -static -fopenmp -std=c++11 -L/usr/local/lib/ -L/usr/local/lib64/ -Wl,--no-as-needed -L${PROJECT_BINARY_DIR}/lib/libloop.a -L/usr/local/lib64/libgomp.a -L/usr/local/lib/libgsl.a -L/usr/local/lib/libgslcblas.a -L/usr/lib64/libblas.a -L/usr/local/lib/libboost_system.a -L/usr/local/lib/libboost_filesystem.a -L/usr/local/lib/boost_timer.a -L/usr/lib64/libpthread_nonshared.a -L/usr/lib64/libc.a -L/usr/lib64/libm.a" )
#SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -std=c++11 -L/usr/local/lib/ -Wl,--no-as-needed -lpthread -lgslcblas -lboost_filesystem -lboost_system -lboost_timer -lgsl -lloop" )
#SET( EXECUTABLE_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/test" )
SET( EXECUTABLE_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/test" )
SET( LIBRARY_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/lib" )
add_library(loop STATIC IMPORTED GLOBAL)
set(C_LIB ${CMAKE_BINARY_DIR}/lib/libloop.a)
set_target_properties(loop
	PROPERTIES
	IMPORTED_LOCATION ${C_LIB}
	)
add_dependencies(loop combinelib)
add_custom_target(combinelib 
        rm -f libloop.a
	COMMAND ar rcT libloop.a libnsppdbstatistics.a libbackbone.a libnsploopclosure.a libnspdataio.a libnspgeometry.a 
		libnspdstl.a libnspkmeans.a libsd.a libnspproteinrep.a libfullsite.a libintegrateddesign.a liballatom.a libnspdesignseq.a #libadn.a
	COMMAND ranlib libloop.a
	DEPENDS nsppdbstatistics backbone nsploopclosure nspdataio nspgeometry nspdstl nspkmeans nspproteinrep sd fullsite 
			integrateddesign allatom nspdesignseq #adn
        WORKING_DIRECTORY ${LIBRARY_OUTPUT_PATH})
add_subdirectory (dstl)
add_subdirectory (dataio)
add_subdirectory (geometry)
add_subdirectory (proteinrep)
add_subdirectory (pdbstatistics)
add_subdirectory (backbone)
add_subdirectory (loopclosure)
add_subdirectory (kmeans)
add_subdirectory (sd)
add_subdirectory (fullsite)
add_subdirectory (designseq)
add_subdirectory (integrateddesign)
#add_subdirectory (addition)
add_subdirectory (allatom)
#add_dependencies (testnestedcontainers combinelib)
#add_subdirectory (app)
