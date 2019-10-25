
project(error_comp)

set(SRCS cases.cpp estimator.cpp utility.cpp projection_factor.cpp pose_local_parameterization.cpp 
	sampson_factor.cpp estimator_depth_factor.cpp projectionOneFrameTwoCamFactor.cpp 
	projectionTwoFrameOneCamFactor.cpp projectionTwoFrameTwoCamFactor.cpp)

add_library(factors ${SRCS})
target_link_libraries(factors ${CERES_LIBRARIES})

add_executable(error_comp main_test.cpp )
target_link_libraries(error_comp factors ${CERES_LIBRARIES})

add_executable(depth_factor_comp main_depth_factor.cpp)
target_link_libraries(depth_factor_comp factors)

add_executable(test_sampson test_sampson.cpp utility.cpp sampson_factor.cpp)
target_link_libraries(test_sampson ${CERES_LIBRARIES})