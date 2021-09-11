
/**
 * \file TriangleQuadrature.h
 * \author Huayi Wei, Chunyu Chen
 *
 * \date 2021.09.08
 *
 * \brief 三角形单元上的积分公式
 *
 */
#ifndef TriangleQuadrature_h
#define TriangleQuadrature_h
#include <array>
#include <vector>

#include "TriangleQuadratureData.h"

namespace OF{
namespace Quadrature{
/**
 * \brief 三角形单元上的积分类
 *
 * \note 这里用 std::array<double, 3> 数组代表积分点.
 *
 */
class TriangleQuadrature
{
public:
  typedef std::array<double, 3> RNode;
public:

  /**
   * \brief 生成积分点和权重
   *
   * \param q 积分的代数精度阶数
   *
   */
  TriangleQuadrature(int q): m_q(q) 
  {
    int N = number_of_quadrature_points();
    m_weight.resize(N);
    m_qpoint.resize(N);
    
    int s = (q*(q-1)+q*(q-1)*(2*q-1)/3)/4;
    for(int i = 0; i < N; i++)
    {
      m_weight[i] = TriangleQuadratureData[i+s][3];
      m_qpoint[i][0] = TriangleQuadratureData[i+s][0];
      m_qpoint[i][1] = TriangleQuadratureData[i+s][1];
      m_qpoint[i][2] = TriangleQuadratureData[i+s][2];
    }
  }

  /**
   * \brief 积分点个数
   */
  int number_of_quadrature_points() 
  {
    return (m_q+1)*m_q/2;
  }

  /**
   *\brief 获取第 i 个积分点
   */
  const RNode & quadrature_point(int i)
  {
    return m_qpoint[i];
  }

  /**
   *\brief 获取第 i 个积分权重
   */
  double quadrature_weight(int i)
  {
    return m_weight[i];
  }

  /**
   *\brief 获取所有积分点
   */
  const std::vector<RNode> & quadrature_points()
  {
    return m_qpoint;
  }

  /**
   *\brief 获取所有积分权重
   */
  const std::vector<double> & quadrature_weights()
  {
    return m_weight;
  }

private:
  int m_q;
  std::vector<double> m_weight; /**< 积分权重 */
  std::vector<RNode>  m_qpoint; /**< 积分点 */
};

};
};
#endif // end of TriangleQuadrature_h
