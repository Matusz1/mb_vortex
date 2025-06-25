project(examples)

add_executable(example_convergence example_convergence.cpp)
target_link_libraries(example_convergence fock::fock)

add_executable(example_energies example_energies.cpp)
target_link_libraries(example_energies fock::fock)

add_executable(example_gpe_bdg example_gpe_bdg.cpp)
target_link_libraries(example_gpe_bdg fock::fock)

add_executable(example_density example_density.cpp)
target_link_libraries(example_density fock::fock)

add_executable(example_particle_drawing example_particle_drawing.cpp)
target_link_libraries(example_particle_drawing fock::fock)

add_executable(example_write_uv example_write_uv.cpp)
target_link_libraries(example_write_uv fock::fock)
