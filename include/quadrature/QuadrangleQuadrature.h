/**
 * \file QuadrangleQuadrature.h 
 * \author Huayi Wei<weihuayi@xtu.edu.cn>
 * \date 2021.09.08
 */
#ifndef QuadrangleQuadrature_h
#define QuadrangleQuadrature_h

#include <array>
#include <vector>

#include "GaussLegendreQuadratureData.h"

namespace OF{
namespace Quadrature{

/**
 * \class QuadrangleQuadrature QuadrangleQuadrature.h
 *
 * \brief \f$[0, 1]^2\f$ 上的积分点和权重
 *
 */
class QuadrangleQuadrature
{
public:
  typedef std::array<double, 2> RNode;
public:

  /**
   * \brief 构造函数
   * \param q 是积分公式的编号
   */
  QuadrangleQuadrature(int q): m_q(q) 
  {
    int NQ = number_of_quadrature_points();
    m_weight.resize(NQ);
    m_qpoint.resize(NQ);
    
    int s = m_q*(m_q+1)/2;
    for(int i = 0; i < m_q+1; i++)
      for(int j = 0; j < m_q+1; j++)
      {
        double w_i = GaussLegendreQuadratureData[s+i][1]/2.0;
        double w_j = GaussLegendreQuadratureData[s+j][1]/2.0;
        m_weight[i] = w_i*w_j;

        int I = i*(m_q+1) + j;
        m_qpoint[I][0] = (GaussLegendreQuadratureData[s+i][0]+1)/2.0; // u
        m_qpoint[I][1] = (GaussLegendreQuadratureData[s+j][0]+1)/2.0; // v
      }
  }

  /**
   * \brief 获取当前积分公式中积分点的个数
   */
  int number_of_quadrature_points() 
  {
    return (m_q+1)*(m_q+1);
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
