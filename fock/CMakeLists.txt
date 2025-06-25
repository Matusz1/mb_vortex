project(mbs)

find_package(GSL REQUIRED)

if (NOT TARGET Eigen)
    find_package(Eigen3 3.4 REQUIRED NO_MODULE)
else()
    add_library(Eigen3::Eigen ALIAS Eigen)
endif()


add_library(fock STATIC
    src/Bogoliubov.cpp
    src/DiskState.cpp
    src/FockCollection.cpp
    src/FockState.cpp
    src/ImportanceTruncation.cpp
    src/Integrator.cpp
    src/MeanFieldSolver.cpp
    src/Numeric.cpp
    src/ParticleDrawing.cpp
    src/Workspace.cpp
)
add_library(fock::fock ALIAS fock)

# Print the main project directory

target_include_directories(fock PUBLIC include)
target_link_libraries(fock
    PRIVATE
        GSL::gsl GSL::gslcblas
    PUBLIC
        Eigen3::Eigen
        Spectra
)
