
/**
 * \file LagrangeFiniteElementSpace3d.h
 * \brief 三维拉格朗日有限元空间
 * \author Chunyu Chen
 * \data 2022/2/13
 */
#ifndef LagrangeFiniteElementSpace3d_h
#define LagrangeFiniteElementSpace3d_h

#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <map>
#include <time.h>

#include <functional>

#include "algebra/linalg.h"
#include "common/Multindex.h"
#include "RArray2d.h"
#include "femcore.h"
#include "barycentric.h"

namespace OF {
namespace FunctionSpace{

/**
 * \brief 第一类棱元空间类
 * \param Mesh 空间的网格类
 * \param AK 代数类
 * \param Quadrature 网格上的积分类, 要有单元积分和边积分功能
 */
template<typename Mesh, typename AK, typename Quadrature>
class LagrangeFiniteElementSpace3d
{
public:
  LagrangeFiniteElementSpace3d(std::shared_ptr<Mesh> & mesh, int p): 
    m_p(p), m_mesh(mesh), m_integrator(p, p+3)
  {
    m_index =   
  }

  void mass_matrix(RArray2d<double> & M)
  {
    mass_matrix(m_integrator->quadrature_points, m_integrator->quadrature_weights, 
        )
  }

private:
  int m_p;
  std::shared_ptr<Mesh> m_mesh;
  std::shared_ptr<Quadrature> m_integrator;
};// end of LagrangeFiniteElementSpace3d

}//end of FunctionSpace
}//end of OF
#endif // end of LagrangeFiniteElementSpace3d_h
