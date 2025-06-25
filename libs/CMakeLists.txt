project(external_libs)

add_library(Spectra INTERFACE)
target_include_directories(Spectra INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/include)

add_library(Eigen INTERFACE)
target_include_directories(Eigen INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Eigen/include)
