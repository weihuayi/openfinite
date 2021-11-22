
configure_file( TestMacro.h.in ${CMAKE_SOURCE_DIR}/test/TestMacro.h)


add_executable(test_TriangleMesh test_TriangleMesh.cpp)
add_test(Mesh test_TriangleMesh)

add_executable(test_ScaledMonomialSpace2d test_ScaledMonomialSpace2d.cpp)
add_test(FemSpace test_ScaledMonomialSpace2d)

add_executable(test_Multindex test_Multindex.cpp)
add_test(Common test_Multindex)

add_executable(test_Permutation test_Permutation.cpp)
add_test(Common test_Permutation)

add_executable(test_IncreaseSeq test_IncreaseSeq.cpp)
add_test(Common test_IncreaseSeq)

add_executable(test_InitTop test_InitTop.cpp)
add_test(InitTop test_InitTop)

add_executable(test_TriangleMeshQuadrature test_TriangleMeshQuadrature.cpp)
add_test(TriangleMeshQuadrature test_TriangleMeshQuadrature)

add_executable(test_array test_array.cpp)
add_test(array test_array)

add_executable(test_Matrix test_Matrix.cpp)
add_test(Algebra test_Matrix)

add_executable(test_Eigen test_Eigen.cpp)
add_test(Algebra test_Eigen)
target_link_libraries(test_Eigen PRIVATE Eigen3::Eigen)

add_executable(test_vector test_vector.cpp)
add_test(Common test_vector)

add_executable(test_map test_map.cpp)
add_test(Common test_map)

add_executable(test_TriangleQuadrature test_TriangleQuadrature.cpp)
add_test(Quadrature test_TriangleQuadrature)

add_executable(test_TetrahedronQuadrature test_TetrahedronQuadrature.cpp)
add_test(Quadrature test_TetrahedronQuadrature)

add_executable(test_GaussLegendreQuadrature test_GaussLegendreQuadrature.cpp)
add_test(Quadrature test_GaussLegendreQuadrature)

add_executable(test_point_in_polygon test_point_in_polygon.cpp)
add_test(PIP test_point_in_polygon)

add_executable(test_Point_2 test_Point_2.cpp)
add_test(Geometry test_Point_2)

if(${OpenBLAS_FOUND})
    add_executable(test_NedecFEM test_NedecFEM.cpp)
    target_link_libraries(test_NedecFEM PRIVATE ${OpenBLAS_LIBRARIES} pthread)
    add_test(FemSpace test_NedecFEM)

    add_executable(test_inv test_inv.cpp)
    target_link_libraries(test_inv PRIVATE ${OpenBLAS_LIBRARIES} pthread)
    add_test(Algebra test_inv)

    add_executable(test_openblas test_openblas.cpp)
    target_link_libraries(test_openblas PRIVATE ${OpenBLAS_LIBRARIES} pthread)
    add_test(Algebra test_openblas)

endif()

if(${VTK_FOUND})
    add_executable(test_HexMesh test_HexMesh.cpp)
    add_test(Mesh test_HexMesh)
    target_link_libraries(test_HexMesh PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS test_HexMesh 
        MODULES ${VTK_LIBRARIES}
    )
    add_executable(test_Dof test_Dof.cpp)
    add_test(FemSpace test_Dof)
    target_link_libraries(test_Dof PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS test_Dof 
        MODULES ${VTK_LIBRARIES}
    )
    add_executable(test_hexmesh_quality test_hexmesh_quality.cpp)
    add_test(MeshQuality test_hexmesh_quality)
    target_link_libraries(test_hexmesh_quality PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS test_hexmesh_quality 
        MODULES ${VTK_LIBRARIES}
    )
    add_executable(test_LagrangeTriangleMesh test_LagrangeTriangleMesh.cpp)
    add_test(Mesh test_LagrangeTriangleMesh)
    target_link_libraries(test_LagrangeTriangleMesh PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS test_LagrangeTriangleMesh 
        MODULES ${VTK_LIBRARIES}
    )
    add_executable(test_LagrangeFESpace2d test_LagrangeFESpace2d.cpp)
    add_test(FemSpace test_LagrangeFESpace2d)
    target_link_libraries(test_LagrangeFESpace2d PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS test_LagrangeFESpace2d 
        MODULES ${VTK_LIBRARIES}
    )
    add_executable(test_BDM test_BDM.cpp)
    add_test(FemSpace test_BDM)
    target_link_libraries(test_BDM PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS test_BDM 
        MODULES ${VTK_LIBRARIES}
    )
endif()


