include(${CMAKE_SOURCE_DIR}/cmake/pr_common.cmake)

set(PROJ_NAME cppbezierfit)
pr_add_library(${PROJ_NAME} STATIC)

pr_add_dependency(${PROJ_NAME} glm_module TARGET PRIVATE)

pr_add_compile_definitions(${PROJ_NAME} -DGLM_ENABLE_EXPERIMENTAL PUBLIC)

pr_init_module(${PROJ_NAME})

pr_finalize(${PROJ_NAME})
