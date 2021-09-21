#ifndef FirstKindNedecFiniteElementSpace2d_h
#define FirstKindNedecFiniteElementSpace2d_h

#include <math.h>
#include <memory>
#include <vector>
#include <array>

#include <functional>

#include "ScaledMonomialSpace2d.h"

namespace OF {
namespace FunctionSpace{

template<typename Mesh>
class FKNDof2d
{
public:
  typedef std::array<int, 2> Dof;

public:
  FKNDof2d(std::shared_ptr<Mesh> mesh, int p): m_mesh(mesh), m_p(p) {}

  int number_of_local_dofs()
  {
    return (m_p+1)*(m_p+3);
  }

  int number_of_edge_local_dofs()
  {
    return m_p+1;
  }

  int number_of_cell_local_dofs()
  {
    return (m_p+1)*m_p;
  }


  int number_of_dofs()
  {
    auto NC = m_mesh->number_of_cells();
    int NE = m_mesh->number_of_edges();
    return NC*(m_p+1)*m_p+NE*(m_p+1);
  }

  /**
   * \brief 第 i 个单元的第 j 个自由度的编号
   */
  int cell_to_dof(int i, int j)
  {
    if(j<3*(m_p+1))
    {
      int idx = j/(m_p+1);
      return (m_p+1)*m_mesh->cell_to_edge(idx)+j-idx*(m_p+1);
    }
    else
    {
      int NE = m_mesh->number_of_edges();
      return NE*(m_p+1)+i*((m_p+1)*m_p)+j-3*(m_p+1);
    }
  }


private:
  int m_p;
  std::shared_ptr<Mesh> m_mesh;
};// end of FKNDof2d

template<typename Mesh, typename AK, typename Quadrature>
class FirstKindNedecFiniteElementSpace2d
{
public:
  typedef FKNDof2d<Mesh> FKNDof;

  typedef typename Mesh::F F;
  typedef typename Mesh::I I;
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;

  typedef typename AK::Matrix Matrix;
  typedef ScaledMonomialSpace2d<Mesh, AK, Quadrature> SMSpace;

public:
  FirstKindNedecFiniteElementSpace2d(std::shared_ptr<Mesh> mesh, int p): 
    m_p(p), m_mesh(mesh), m_integrator(mesh, p+3)
  {
    m_smspace = std::make_shared<SMSpace>(mesh, p); 
    m_dof = std::make_shared<FKNDof>(mesh, p); /**< \todo 为什么放在初始化里不行 */
  }

  void basis_coefficients(int idx, Matrix & coeff)
  {
    auto ldof = m_dof->number_of_local_dofs();
    coeff.reshape(ldof, ldof);
    auto cell = m_mesh->cell(idx);
    auto edge = m_mesh->cell_to_edge(idx); /**< 获取三条边 */

    auto eldof = m_p+1;
    auto cldof = (m_p+1)*m_p;
    int smdof = (m_p+1)*(m_p+2)/2;

    for(int i = 0; i < 3; i++)
    {
      Matrix mat;/**< 边上的积分 */
      int e = edge[i];

      /**
       * 定义被积函数
       */
      std::function<void(const Node&, Matrix&)> fm = [this, &idx, &e, &ldof, &cldof, &eldof, &smdof]
        (const Node & p, Matrix & mat)->void
      {
        auto ldof = m_dof->number_of_local_dofs();
        if(mat.number_of_rows()!=m_p+1 || mat.number_of_columns()!=ldof)
        {
          mat.reshape(m_p+1, ldof);
        }

        std::vector<F> mval, fval; /**< mval: 单元上的基函数在 p 处的值, fval 边上的基函数在 p 处的值 */
        m_smspace->basis(idx, p, mval, m_p>1?m_p:1);
        m_smspace->edge_basis(e, p, fval);

        auto t = m_mesh->edge_tangent(e);

        for(int j = 0; j < eldof; j++)
        {
          for(int k = 0; k < smdof; k++)
          {
            auto val = fval[j]*mval[k];
            mat[j][k] = t[0]*val;
            mat[j][k+smdof] = t[1]*val;
          }
          for(int k = smdof*2; k < ldof; k++)
          {
            mat[j][k] = (t[0]*mval[2]-t[1]*mval[1])*fval[j]*mval[k-smdof*2+cldof/2];
          }
        }
      }; /**< 函数定义完成 */
      m_integrator.edge_integral(e, fm, mat); /**< 将函数积分, 得到的矩阵 mat 是 coeff 的一部分 */
      for(int j = 0; j < eldof; j++)
      {
        for(int k = 0; k < ldof; k++)
        {
          coeff[j+i*eldof][k] = mat[j][k];
        }
      }
    }
    Matrix mass;
    m_smspace->mass_matrix(idx, mass);

    /**
     * \note 生成 x*m_k 的编号
     * \todo 用更优雅的方式实现
     */
    int xidx[cldof];
    int x = 0, index = 0;
    for(int i = 0; i < m_p+1; i++)
    {
      for(int j = 0; j < i; j++)
      {
        xidx[index++] = x++; 
      }
      x++;
    }

    for(int i = 0; i < cldof/2; i++)
    {
      for(int j = 0; j < smdof; j++)
      {
        coeff[i+3*eldof][j] = mass[i][j]; 
        coeff[i+3*eldof+cldof/2][j+smdof] = mass[i][j]; 
      }
      for(int j = 0; j < eldof; j++)
      {
        coeff[i+3*eldof][j+smdof*2] = mass[xidx[i]+1][j+cldof/2]; 
        coeff[i+3*eldof+cldof/2][j+smdof*2] = -mass[xidx[i]][j+cldof/2]; 
      }
    }
  }

private:
  int m_p;
  Quadrature m_integrator;
  std::shared_ptr<FKNDof> m_dof;
  std::shared_ptr<Mesh> m_mesh;
  std::shared_ptr<SMSpace> m_smspace;
};// end of FirstKindNedecFiniteElementSpace2d

}//end of FunctionSpace
}//end of OF
#endif // end of FirstKindNedecFiniteElementSpace2d_h
