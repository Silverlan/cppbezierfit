include(${CMAKE_SOURCE_DIR}/cmake/pr_common.cmake)

set(PROJ_NAME cppbezierfit)
pr_add_library(${PROJ_NAME} STATIC)

pr_add_external_dependency(${PROJ_NAME} glm HEADER_ONLY PRIVATE)

pr_init_module(${PROJ_NAME})

pr_finalize(${PROJ_NAME})
