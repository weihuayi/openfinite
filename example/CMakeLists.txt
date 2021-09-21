
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
endif()

if(${CGAL_FOUND} AND ${VTK_FOUND})
    add_executable(pmesh_tetmesh_generation pmesh_tetmesh_generation.cpp)
    target_link_libraries(pmesh_tetmesh_generation PRIVATE ${CGAL_LIBRARIES}
        ${VTK_LIBRARIES} ${METIS_LIBRARY})
    vtk_module_autoinit(
        TARGETS pmesh_tetmesh_generation 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(pmesh_trimesh_generation pmesh_trimesh_generation.cpp)
    target_link_libraries(pmesh_trimesh_generation PRIVATE ${CGAL_LIBRARIES}
        ${VTK_LIBRARIES} ${METIS_LIBRARY})
    vtk_module_autoinit(
        TARGETS pmesh_trimesh_generation 
        MODULES ${VTK_LIBRARIES}
    )
endif()

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


endif()
