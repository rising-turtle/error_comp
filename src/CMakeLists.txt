

set(SRCS cases.cpp estimator.cpp utility.cpp projection_factor.cpp pose_local_parameterization.cpp 
	sampson_factor.cpp main_test.cpp)

add_executable(error_comp ${SRCS})
target_link_libraries(error_comp ${CERES_LIBRARIES})

add_executable(test_sampson test_sampson.cpp utility.cpp sampson_factor.cpp)
target_link_libraries(test_sampson ${CERES_LIBRARIES})