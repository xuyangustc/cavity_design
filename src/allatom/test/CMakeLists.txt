SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
SET(BUILD_SHARED_LIBS OFF)
SET(CMAKE_EXE_LINKER_FLAGS "-static")

link_directories(/usr/local/lib)
link_directories(/usr/local/lib64)
link_directories(/usr/lib)
link_directories(/usr/lib64)

link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable (Initial_Chain_Gen cav_initial_chain_gen.cpp)
add_executable (Loop_replace cav_loop_replace_round_1.cpp)
add_executable (Loop_replace_2 cav_loop_replace_round_2.cpp)
add_executable (SDrun_far cav_sd_fix_att_rep.cpp)
add_executable (SDrun cav_sd.cpp)

target_link_libraries (Initial_Chain_Gen loop gomp gsl gslcblas blas boost_system boost_filesystem boost_timer pthread_nonshared m)
target_link_libraries (Loop_replace loop gomp gsl gslcblas blas boost_system boost_filesystem boost_timer pthread_nonshared m)
target_link_libraries (Loop_replace_2 loop gomp gsl gslcblas blas boost_system boost_filesystem boost_timer pthread_nonshared m)
target_link_libraries (SDrun_far loop gomp gsl gslcblas blas boost_system boost_filesystem boost_timer pthread_nonshared m)
target_link_libraries (SDrun loop gomp gsl gslcblas blas boost_system boost_filesystem boost_timer pthread_nonshared m)
