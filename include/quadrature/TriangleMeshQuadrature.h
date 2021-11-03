/**
 * \file TriangleMeshQuadrature.h
 * \author 陈春雨
 * \date 2021/9/11
 *
 * \brief 三角形网格上积分
 */

#ifndef TriangleMeshQuadrature_h
#define TriangleMeshQuadrature_h

#include <math.h>
#include <memory>
#include <functional>

#include "TriangleQuadrature.h"
#include "GaussLegendreQuadrature.h"

namespace OF{
namespace Quadrature{

/**
 * \brief 三角形网格上的积分类
 * \param Mesh 是任何一种三角形网格
 */
template<typename Mesh>
class TriangleMeshQuadrature
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
  TriangleMeshQuadrature(std::shared_ptr<Mesh> mesh, int p): 
    m_mesh(mesh), m_integrator(p), m_edgeintegrator(p)
  {
    mesh->cell_measure(m_cellmeasure);
    mesh->cell_barycenter(m_cellbarycenter);
    mesh->edge_measure(m_edgemeasure);
    mesh->edge_barycenter(m_edgebarycenter);
  }

  /**
   * \brief 积分函数, 对于被给函数在被给单元上积分, 返回积分值
   * \param i 积分单元
   * \param f 被积函数, 要求是 Point -> R 的函数
   */
  F integral(int i, std::function<F(const Node &, I)> f)
  {
    auto & cellarea = m_cellmeasure[i]; 
    auto & cellbc = m_cellbarycenter[i];
    auto & cell = m_mesh->cell(i);

    auto & p0 = m_mesh->node(cell[0]);
    auto & p1 = m_mesh->node(cell[1]);
    auto & p2 = m_mesh->node(cell[2]);

    int N = m_integrator.number_of_quadrature_points();
    double val = 0.0;
    for(int j = 0; j < N; j++)
    {
      auto w = m_integrator.quadrature_weight(j);/**< 积分点的权重 */
      auto & qpts = m_integrator.quadrature_point(j); /**< 获取积分点重心坐标 */
      auto P = qpts[0]*p0 + qpts[1]*p1 + qpts[2]*p2; 
      val += f(P, i)*w;
    }
    val *= cellarea;
    return val;
  }

  /**
   * \brief 积分函数, 对于被给函数在被给网格上积分, 返回积分值
   * \param i 积分单元
   * \param f 被积函数, 要求是 Point, I -> R 的函数
   */
  F integral(std::function<F(const Node &, I)> f)
  {
    double val = 0.0;
    auto NC = m_mesh->number_of_cells();
    for(int i = 0; i < NC; i++)
      val += integral(i, f);
    return val;
  }

  /**
   * \brief 对矩阵函数一个单元上积分
   * \param i 积分单元
   * \param f 被积函数, 要求是 Point -> R 的函数
   * \param Tensor 模版参数, 可以是任意有 "*(double)", "*=(double)", "+=(Tensor)" 
   * 运算的数据结构
   */
  template<typename Tensor>
  void integral(int i, std::function<void(const Node&, Tensor&)> & f, Tensor & mat)
  {
    auto & cellarea = m_cellmeasure[i]; 
    auto & cellbc = m_cellbarycenter[i];
    auto & cell = m_mesh->cell(i);

    auto & p0 = m_mesh->node(cell[0]);
    auto & p1 = m_mesh->node(cell[1]);
    auto & p2 = m_mesh->node(cell[2]);

    int N = m_integrator.number_of_quadrature_points();

    auto w = m_integrator.quadrature_weight(0);/**< 积分点的权重 */
    auto & qpts = m_integrator.quadrature_point(0); /**< 获取积分点重心坐标 */
    auto P = qpts[0]*p0 + qpts[1]*p1 + qpts[2]*p2; 
    f(P, mat);
    mat *= w;
    for(int j = 1; j < N; j++)
    {
      auto w = m_integrator.quadrature_weight(j);/**< 积分点的权重 */
      auto & qpts = m_integrator.quadrature_point(j); /**< 获取积分点重心坐标 */
      auto P = qpts[0]*p0 + qpts[1]*p1 + qpts[2]*p2; 
      Tensor tmp;
      f(P, tmp);
      mat += w*tmp;
    }
    mat *= cellarea;
  }

  /**
   * \brief 对矩阵函数一条边上积分
   * \param i 积分单元
   * \param f 被积函数, 要求是 Point -> R 的函数
   * \param Tensor 模版参数, 可以是任意有 "*(double)", "*=(double)", "+=(Tensor)" 
   * 运算的数据结构
   */
  template<typename Tensor>
  void edge_integral(int i, std::function<void(const Node&, Tensor&)> & f, Tensor & mat)
  {
    auto & edgearea = m_edgemeasure[i]; 
    auto & edgebc = m_edgebarycenter[i];
    auto & edge = m_mesh->edge(i);

    auto & p0 = m_mesh->node(edge[0]);
    auto & p1 = m_mesh->node(edge[1]);

    int N = m_edgeintegrator.number_of_quadrature_points();
    auto w = m_edgeintegrator.quadrature_weight(0);/**< 积分点的权重 */
    auto & qpts = m_edgeintegrator.quadrature_point(0); /**< 获取积分点重心坐标 */
    auto P = qpts[0]*p0 + qpts[1]*p1;
    f(P, mat);
    mat *= w;
    for(int j = 1; j < N; j++)
    {
      auto w = m_edgeintegrator.quadrature_weight(j);/**< 积分点的权重 */
      auto & qpts = m_edgeintegrator.quadrature_point(j); /**< 获取积分点重心坐标 */
      auto P = qpts[0]*p0 + qpts[1]*p1;
      Tensor tmp;
      f(P, tmp);
      mat += w*tmp;
    }
    mat *= edgearea;
  }

private:
  std::shared_ptr<Mesh> m_mesh; /**< 积分区域的网格 */
  std::vector<F> m_cellmeasure; /**< 网格每个单元的面积 */
  std::vector<Node> m_cellbarycenter; /**< 网格中每个单元的重心 */
  TriangleQuadrature m_integrator; /**< 单元上的积分算法 */
  std::vector<F> m_edgemeasure; /**< 网格每个单元的面积 */
  std::vector<Node> m_edgebarycenter; /**< 网格中每个单元的重心 */
  GaussLegendreQuadrature m_edgeintegrator; /**< 单元上的积分算法 */
};

}
}


#endif // end of TriangleMeshQuadrature_h
