
configure_file( TestMacro.h.in ${CMAKE_SOURCE_DIR}/test/TestMacro.h)


add_executable(test_TriangleMesh test_TriangleMesh.cpp)
add_test(Mesh test_TriangleMesh)

add_executable(test_ScaledMonomialSpace2d test_ScaledMonomialSpace2d.cpp)
add_test(FemSpace test_ScaledMonomialSpace2d)

add_executable(test_InitTop test_InitTop.cpp)
add_test(InitTop test_InitTop)

add_executable(test_TriangleMeshQuadrature test_TriangleMeshQuadrature.cpp)
add_test(TriangleMeshQuadrature test_TriangleMeshQuadrature)

add_executable(test_array test_array.cpp)
add_test(array test_array)

add_executable(test_Matrix test_Matrix.cpp)
add_test(Algebra test_Matrix)

add_executable(test_vector test_vector.cpp)
add_test(Common test_vector)

add_executable(test_TriangleQuadrature test_TriangleQuadrature.cpp)
add_test(Quadrature test_TriangleQuadrature)

add_executable(test_TetrahedronQuadrature test_TetrahedronQuadrature.cpp)
add_test(Quadrature test_TetrahedronQuadrature)

add_executable(test_GaussLegendreQuadrature test_GaussLegendreQuadrature.cpp)
add_test(Quadrature test_GaussLegendreQuadrature)

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

endif()


