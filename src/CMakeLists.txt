

set(SRCS cases.cpp estimator.cpp utility.cpp projection_factor.cpp pose_local_parameterization.cpp 
	main_test.cpp)

add_executable(error_comp ${SRCS})
target_link_libraries(error_comp ${CERES_LIBRARIES})