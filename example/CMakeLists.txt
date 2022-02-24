
if(${VTK_FOUND})

    add_executable(pmesh_hexmesh_generation pmesh_hexmesh_generation.cpp)
    target_link_libraries(pmesh_hexmesh_generation PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS pmesh_hexmesh_generation 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(pmesh_generation pmesh_generation.cpp)
    target_link_libraries(pmesh_generation PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS pmesh_generation 
        MODULES ${VTK_LIBRARIES}
    )

#    add_executable(pmesh_tetmesh_generation pmesh_tetmesh_generation.cpp)
#    target_link_libraries(pmesh_tetmesh_generation PRIVATE ${CGAL_LIBRARIES}
#        ${VTK_LIBRARIES} ${METIS_LIBRARY})
#    vtk_module_autoinit(
#        TARGETS pmesh_tetmesh_generation 
#        MODULES ${VTK_LIBRARIES}
#    )
#
#    add_executable(pmesh_trimesh_generation pmesh_trimesh_generation.cpp)
#    target_link_libraries(pmesh_trimesh_generation PRIVATE ${CGAL_LIBRARIES}
#        ${VTK_LIBRARIES} ${METIS_LIBRARY})
#    vtk_module_autoinit(
#        TARGETS pmesh_trimesh_generation 
#        MODULES ${VTK_LIBRARIES}
#    )
#
    add_executable(functional_programming_fem_example
        functional_programming_fem_example.cpp)
    target_link_libraries(functional_programming_fem_example PRIVATE ${CGAL_LIBRARIES}
        ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS functional_programming_fem_example 
        MODULES ${VTK_LIBRARIES}
    )

endif()

#if(${CGAL_FOUND} AND ${VTK_FOUND})
#endif()

if(${VTK_FOUND} AND ${MPI_FOUND})

    add_executable(parallel_tetmesh_optimization parallel_tetmesh_optimization.cpp)
    target_link_libraries(parallel_tetmesh_optimization PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_tetmesh_optimization 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_quadmesh_optimization parallel_quadmesh_optimization.cpp)
    target_link_libraries(parallel_quadmesh_optimization PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_quadmesh_optimization 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_trimesh_optimization parallel_trimesh_optimization.cpp)
    target_link_libraries(parallel_trimesh_optimization PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_trimesh_optimization 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_hexmesh_optimization parallel_hexmesh_optimization.cpp)
    target_link_libraries(parallel_hexmesh_optimization PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_hexmesh_optimization 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_tetmesh_optimization0 parallel_tetmesh_optimization0.cpp)
    target_link_libraries(parallel_tetmesh_optimization0 PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_tetmesh_optimization0 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_tetmesh_optimization1 parallel_tetmesh_optimization1.cpp)
    target_link_libraries(parallel_tetmesh_optimization1 PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_tetmesh_optimization1 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_quadmesh_optimization0 parallel_quadmesh_optimization0.cpp)
    target_link_libraries(parallel_quadmesh_optimization0 PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_quadmesh_optimization0 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_trimesh_optimization0 parallel_trimesh_optimization0.cpp)
    target_link_libraries(parallel_trimesh_optimization0 PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_trimesh_optimization0 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_hexmesh_optimization0 parallel_hexmesh_optimization0.cpp)
    target_link_libraries(parallel_hexmesh_optimization0 PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_hexmesh_optimization0 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_hexmesh_optimization1 parallel_hexmesh_optimization1.cpp)
    target_link_libraries(parallel_hexmesh_optimization1 PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_hexmesh_optimization1
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_mesh_coloring parallel_mesh_coloring.cpp)
    target_link_libraries(parallel_mesh_coloring PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_mesh_coloring 
        MODULES ${VTK_LIBRARIES}
    )
    add_executable(Poisson Poisson.cpp )
    target_link_libraries(Poisson PRIVATE ${VTK_LIBRARIES} Eigen3::Eigen)
    vtk_module_autoinit(
        TARGETS Poisson
        MODULES ${VTK_LIBRARIES}
    )
    add_executable(PoissonWithMix PoissonWithMix.cpp )
    target_link_libraries(PoissonWithMix PRIVATE ${VTK_LIBRARIES} Eigen3::Eigen)
    vtk_module_autoinit(
        TARGETS PoissonWithMix
        MODULES ${VTK_LIBRARIES}
    )

endif()
