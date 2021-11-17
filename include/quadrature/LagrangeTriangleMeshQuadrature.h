
/**
 * \file LagrangeTriangleMeshQuadrature.h
 * \author 陈春雨
 * \date 2021/9/11
 *
 * \brief 三角形网格上积分
 */

#ifndef LagrangeTriangleMeshQuadrature_h
#define LagrangeTriangleMeshQuadrature_h

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
class LagrangeTriangleMeshQuadrature
{
public:
  typedef typename Mesh::I I;
  typedef typename Mesh::F F;

public:
  /**
   * \brief 构造函数
   * \param mesh 积分区域, 是一个智能指针
   * \param p 积分的代数精度
   */
  LagrangeTriangleMeshQuadrature(std::shared_ptr<Mesh> mesh, int p): 
    m_mesh(mesh), m_integrator(p), m_edgeintegrator(p) {}

  /**
   * \brief 积分函数, 对于被给函数在被给单元上积分, 返回积分值
   * \param i 积分单元
   * \param f 被积函数, 要求是 Point -> R 的函数
   */
  template<typename Container>
  F integral(int i, std::function<F(const Container &)> f)
  {
    int N = m_integrator.number_of_quadrature_points();
    double val = 0.0;
    for(int j = 0; j < N; j++)
    {
      auto w = m_integrator.quadrature_weight(j);/**< 积分点的权重 */
      auto & qpts = m_integrator.quadrature_point(j); /**< 获取积分点重心坐标 */

      auto J = m_mesh->det_jacobi(i, qpts);
      val += f(qpts)*w*J;
    }
    val /= 2.0;
    return val;
  }

  /**
   * \brief 积分函数, 对于被给函数在被给网格上积分, 返回积分值
   * \param i 积分单元
   * \param f 被积函数, 要求是 (Container, I) -> R 的函数
   */
  template<typename Container>
  F integral(std::function<F(const Container &, I)> f)
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
  template<typename Tensor, typename Container>
  void integral(int i, std::function<void(const Container&, Tensor&)> & f, Tensor & mat)
  {
    int N = m_integrator.number_of_quadrature_points();
    auto w = m_integrator.quadrature_weight(0);/**< 积分点的权重 */
    auto & qpts = m_integrator.quadrature_point(0); /**< 获取积分点重心坐标 */
    f(qpts, mat);
    auto J = m_mesh->det_jacobi(i, qpts);
    mat *= w*J;
    for(int j = 1; j < N; j++)
    {
      auto w = m_integrator.quadrature_weight(j);/**< 积分点的权重 */
      auto & qpts = m_integrator.quadrature_point(j); /**< 获取积分点重心坐标 */
      Tensor tmp;
      f(qpts, tmp);
      J = m_mesh->det_jacobi(i, qpts);
      mat += w*tmp*J;
    }
    mat = mat*0.5;
  }

  /**
   * \brief 对矩阵函数一条边上积分
   * \param i 积分单元
   * \param f 被积函数, 要求是 Point -> R 的函数
   * \param Tensor 模版参数, 可以是任意有 "*(double)", "*=(double)", "+=(Tensor)" 
   * 运算的数据结构
   * \todo TODO 积分需要测试正确性.
   */
  template<typename Tensor, typename Container>
  void edge_integral(int i, std::function<void(const Container&, Tensor&)> & f, Tensor & mat)
  {
    int N = m_edgeintegrator.number_of_quadrature_points();
    auto w = m_edgeintegrator.quadrature_weight(0);/**< 积分点的权重 */
    auto & qpts = m_edgeintegrator.quadrature_point(0); /**< 获取积分点重心坐标 */
    f(qpts, mat);
    mat *= w;
    for(int j = 1; j < N; j++)
    {
      auto w = m_edgeintegrator.quadrature_weight(j);/**< 积分点的权重 */
      auto & qpts = m_edgeintegrator.quadrature_point(j); /**< 获取积分点重心坐标 */
      Tensor tmp;
      f(qpts, tmp);
      mat += w*tmp;
    }
  }

private:
  std::shared_ptr<Mesh> m_mesh; /**< 积分区域的网格 */
  TriangleQuadrature m_integrator; /**< 单元上的积分算法 */
  GaussLegendreQuadrature m_edgeintegrator; /**< 单元上的积分算法 */
};

}
}


#endif // end of LagrangeTriangleMeshQuadrature_h
