/**
 * \file HexahedronQuadrature.h 
 * \author Huayi Wei<weihuayi@xtu.edu.cn>
 * \date 2021.09.08
 */
#ifndef HexahedronQuadrature_h
#define HexahedronQuadrature_h

#include <array>
#include <vector>

#include "GaussLegendreQuadratureData.h"

namespace OF{
namespace Quadrature{

/**
 * \class HexahedronQuadrature HexahedronQuadrature.h
 *
 * \brief \f$[0, 1]^3\f$ 上的积分点和权重
 *
 * 注意这里给的是积分点处的形函数值
 *
 * \f[
 *    (1-u)*(1-v)*(1-w), u*(1-v)*(1-w), u*v*(1-w), (1-u)*v*(1-w),
 *    (1-u)*(1-v)*w, u*(1-v)*w, u*v*w, (1-u)*v*w 
 * \f]
 *
 * 其中 \f$(u, v, w) \in [0, 1]^3\f$
 */
class HexahedronQuadrature
{
public:
  typedef std::array<double, 3> RNode;

public:

  /**
   * \brief 构造函数
   * \param q 
   */
  HexahedronQuadrature(int q): m_q(q) 
  {
    int NQ = number_of_quadrature_points();
    m_weight.resize(NQ);
    m_qpoint.resize(NQ);
    
    int s = m_q*(m_q+1)/2;
    for(int i = 0; i < m_q+1; i++)
      for(int j = 0; j < m_q+1; j++)
        for(int k = 0; k < m_q+1; k++) 
      {
        double w_i = GaussLegendreQuadratureData[s+i][1]/2.0;
        double w_j = GaussLegendreQuadratureData[s+j][1]/2.0;
        double w_k = GaussLegendreQuadratureData[s+k][1]/2.0;
        m_weight[i] = w_i*w_j*w_k;

        int I = i*(m_q + 1)*(m_q + 1) + j*(m_q + 1) + k;
        m_qpoint[I][0] = (GaussLegendreQuadratureData[s+i][0]+1)/2.0;//u
        m_qpoint[I][1] = (GaussLegendreQuadratureData[s+j][0]+1)/2.0;//v
        m_qpoint[I][2] = (GaussLegendreQuadratureData[s+k][0]+1)/2.0;//w

      }
  }

  /**
   * \brief 获取当前积分公式中积分点的个数
   */
  int number_of_quadrature_points() 
  {
    return (m_q+1)*(m_q+1)*(m_q+1);
  } 

  /**
   * \brief 获取当前积分公式中第 i 个参考积分点
   */
  const RNode & quadrature_point(int i)
  {
    return m_qpoint[i];
  }

  const double & quadrature_weight(int i)
  {
    return m_weight[i];
  }

  const std::vector<RNode> & quadrature_points()
  {
    return m_qpoint;
  }

  const std::vector<double> & quadrature_weights()
  {
    return m_weight;
  }

private:
  int m_q;
  std::vector<double> m_weight;
  std::vector<RNode>  m_qpoint;
};

} // end of namespace Quadrature
} // end of namespace TOPT
#endif // end of QuadrangleQuadrature_h
