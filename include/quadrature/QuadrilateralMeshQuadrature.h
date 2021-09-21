
 /**
 * \file QuadrilateralMeshQuadrature.h
 * \author 陈春雨
 * \date 2021/9/11
 *
 * \brief 三角形网格上积分
 */

#ifndef QuadrilateralMeshQuadrature_h
#define QuadrilateralMeshQuadrature_h

#include <math.h>
#include <memory>
#include <functional>


namespace OF{
namespace Quadrature{

/**
 * \brief 三角形网格上的积分类
 * \param Mesh 是任何一种三角形网格
 */
template<typename Mesh>
class QuadrilateralMeshQuadrature
{
public:
  typedef typename Mesh::I I;
  typedef typename Mesh::F F;
  typedef typename Mesh::Node Node;

public:
  /**
   * \brief 构造函数
   * \param mesh 积分区域, 是一个智能指针
   * \param p 积分的代数精度
   */
  QuadrilateralMeshQuadrature(std::shared_ptr<Mesh> mesh, int p): m_mesh(mesh)
  {
    mesh->cell_measure(m_cellmeasure);
    mesh->cell_barycenter(m_cellbarycenter);
  }

private:
  std::shared_ptr<Mesh> m_mesh; /**< 积分区域的网格 */
  std::vector<F> m_cellmeasure; /**< 网格每个单元的面积 */
  std::vector<Node> m_cellbarycenter; /**< 网格中每个单元的重心 */
};

}
}


#endif // end of QuadrilateralMeshQuadrature_h
