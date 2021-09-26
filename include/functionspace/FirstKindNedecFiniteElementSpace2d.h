/**
 * \file FirstKindNedecFiniteElementSpace2d.h
 * \brief 二维第一类棱元空间
 * \author 魏华祎, 陈春雨
 * \data 2021/09/15
 */
#ifndef FirstKindNedecFiniteElementSpace2d_h
#define FirstKindNedecFiniteElementSpace2d_h

#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <map>

#include <functional>

#include "algebra/linalg.h"
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
      int eidx = m_mesh->cell_to_edge(i)[idx];
      return (m_p+1)*eidx+j-idx*(m_p+1);
    }
    else
    {
      int NE = m_mesh->number_of_edges();
      return NE*(m_p+1)+i*((m_p+1)*m_p)+j-3*(m_p+1);
    }
  }

  /**
   * \brief 第 i 个单元的自由度的编号
   */
  std::vector<int> cell_to_dof(int i)
  {
    auto ldof = number_of_local_dofs();
    auto cldof = number_of_cell_local_dofs();
    auto eldof = number_of_edge_local_dofs();
    std::vector<int> celldof(ldof);
    int NE = m_mesh->number_of_edges();
    auto e = m_mesh->cell_to_edge(i);
    for(int j = 0; j < eldof; j++)
    {
      celldof[j] = e[0]*eldof+j;
      celldof[j+eldof] = e[1]*eldof+j;
      celldof[j+eldof*2] = e[2]*eldof+j;
    }
    for(int j = 0; j < cldof; j++)
    {
      celldof[j+3*eldof] = NE*eldof+i*cldof+j;
    }
    return celldof;
  }

private:
  int m_p;
  std::shared_ptr<Mesh> m_mesh;
};// end of FKNDof2d


/**
 * \brief 第一类棱元空间类
 * \param Mesh 空间的网格类
 * \param AK 代数类
 * \param Quadrature 网格上的积分类, 要有单元积分和边积分功能
 */
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
  typedef typename std::vector<std::map<I, F> > BSMatrix;
  typedef ScaledMonomialSpace2d<Mesh, AK, Quadrature> SMSpace;

public:

  /**
   * \brief 有限元空间构造函数
   * \param mesh 空间的网格
   * \param p 空间次数
   * \todo m_dof 和 m_space 为什么放在初始化里不行
   * \todo 积分子有两个
   */
  FirstKindNedecFiniteElementSpace2d(std::shared_ptr<Mesh> mesh, int p): 
    m_p(p), m_mesh(mesh), m_integrator(mesh, p+3)
  {
    m_smspace = std::make_shared<SMSpace>(mesh, p); 
    m_dof = std::make_shared<FKNDof>(mesh, p);
    build_coeff();
  }

  /**
   * \brief 构建基函数的系数 
   */
  void build_coeff()
  {
    auto NC = m_mesh->number_of_cells();
    m_coeff.resize(NC);
    for(int i = 0; i < NC; i++)
    {
      basis_coefficients(i, m_coeff[i]);
    }
  }

  /**
   * \brief 基函数的系数, 每个单元有一组可以计算的函数 Phi, 求有限元基函数在 Phi 上的系数
   * \param idx 单元编号
   * \param coeff 有限元基函数在 Phi 中的系数, 每一列表示一个基函数.
   */
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

      /** 定义被积函数 */
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

    /** 求逆 */
    OF::AlgebraAlgrithom::matInv(coeff);
  }

  /**
   * \brief 一个单元上基函数在一个点上的值
   * \param cidx 单元编号
   * \param point 节点的笛卡尔坐标
   * \param phi 返回值.
   */
  void basis(int cidx, const Node & point, std::vector<Vector> & phi)
  {
    /** 获得基函数系数 */
    const auto & coeff = m_coeff[cidx];
    std::vector<double> smphi;
    m_smspace->basis(cidx, point, smphi, m_p+1);

    auto ldof = m_dof->number_of_local_dofs();
    auto cldof = m_dof->number_of_cell_local_dofs();
    auto eldof = m_dof->number_of_edge_local_dofs();
    auto smdof = (m_p+2)*(m_p+1)/2;

    phi.resize(ldof);
    for(int i = 0; i < ldof; i++)
    {
      phi[i] = Vector(0, 0); /**< 初始化测试 */
      for(int j = 0; j < smdof; j++)
      {
        phi[i][0] += coeff[j][i]*smphi[j];
        phi[i][1] += coeff[j+smdof][i]*smphi[j];
      }
      for(int j = 0; j < eldof; j++)
      {
        phi[i][0] += coeff[j+smdof*2][i]*smphi[smdof+j+1];
        phi[i][1] += -coeff[j+smdof*2][i]*smphi[smdof+j];
      }
    }
  }

  /**
   * \brief 一个单元上基函数的旋度在一个点上的值
   * \param cidx 单元编号
   * \param point 节点的笛卡尔坐标
   * \param phi 返回值.
   */
  void curl_basis(int cidx, const Node & point, std::vector<double> & phi)
  {
    /** 获得基函数系数 */
    const auto & coeff = m_coeff[cidx];
    std::vector<Vector> smgphi;
    m_smspace->grad_basis(cidx, point, smgphi, m_p+1);

    auto ldof = m_dof->number_of_local_dofs();
    auto cldof = m_dof->number_of_cell_local_dofs();
    auto eldof = m_dof->number_of_edge_local_dofs();
    auto smdof = (m_p+2)*(m_p+1)/2;

    phi.resize(ldof);
    for(int i = 0; i < ldof; i++)
    {
      phi[i] = 0.0; /**< 初始化测试 */
      for(int j = 0; j < smdof; j++)
      {
        phi[i] += -coeff[j][i]*smgphi[j][1];
        phi[i] += coeff[j+smdof][i]*smgphi[j][0];
      }
      for(int j = 0; j < eldof; j++)
      {
        phi[i] += -coeff[j+smdof*2][i]*(smgphi[smdof+j+1][1]+smgphi[smdof+j][0]);
      }
    }
  }

  /**
   * \brief 单元上的基函数的基函数在一点处的值
   * \param cidx 单元编号
   * \param point 节点的笛卡尔坐标
   * \param smvphi 返回值
   */
  void _basis(int cidx, const Node & point, std::vector<Vector> & smvphi)
  {
    std::vector<double> smphi;
    m_smspace->basis(cidx, point, smphi, m_p+1);

    auto ldof = m_dof->number_of_local_dofs();
    auto cldof = m_dof->number_of_cell_local_dofs();
    auto eldof = m_dof->number_of_edge_local_dofs();
    auto smdof = (m_p+2)*(m_p+1)/2;

    smvphi.resize(ldof);
    for(int i = 0; i < smdof; i++)
    {
      smvphi[i][0] = smphi[i]; 
      smvphi[i][1] = 0.0; 
      smvphi[i][0] = smphi[i]; 
      smvphi[i][1] = 0.0; 
    }
    for(int i = 0; i < eldof; i++)
    {
      smvphi[i+smdof*2][0] = smphi[smdof+i+1];
      smvphi[i+smdof*2][1] = -smphi[smdof+i];
    }
  }

  /**
   * \brief 单元上的基函数的质量矩阵
   * \param cidx 单元编号
   * \param mass 返回的质量矩阵
   */
  void _mass_matrix(int cidx, Matrix & mass)
  {
    Matrix smMass;
    m_smspace->mass_matrix(cidx, smMass, m_p+1);

    auto ldof = m_dof->number_of_local_dofs();
    mass.reshape(ldof, ldof);
    int smdof = (m_p+1)*(m_p+2)/2;
    auto eldof = m_dof->number_of_edge_local_dofs();
    for(int i = 0; i < smdof; i++)
    {
      for(int j = 0; j < smdof; j++)
      {
        mass[i][j] = smMass[i][j];
        mass[i+smdof][j+smdof] = smMass[i][j];
      }
      for(int j = 0; j < eldof; j++)
      {
        mass[i][j+smdof*2] = smMass[i][smdof+j+1];
        mass[i+smdof][j+smdof*2] = -smMass[i][smdof+j];
      }
    }
    for(int i = smdof*2; i < ldof; i++)
    {
      for(int j = 0; j < smdof*2; j++)
      {
        mass[i][j] = mass[j][i];
      }
      for(int j = 0; j < eldof; j++)
      {
        mass[i][j+smdof*2] = smMass[-smdof+i+1][smdof+j+1]+smMass[-smdof+i][smdof+j];
      }
    }
  }

  /**
   * \brief 单元上的质量矩阵
   * \param cidx 单元编号
   * \param mass 返回的质量矩阵
  void mass_matrix(int cidx, Matrix & mass)
  {
    auto & coeff = m_coeff[cidx];
    Matrix _mass;
    _mass_matrix(cidx, _mass);
    mass = coeff.transpose_multiply(_mass*coeff);
  }
   */
  void mass_matrix(int i, Matrix & mass)
  {
    auto ldof = m_dof->number_of_local_dofs();
    std::function<void(const Node&, Matrix&)> phi_jk = [this, &i]
      (const Node & point, Matrix & m)->void
    {
      std::vector<Vector> val;
      basis(i, point, val);
      int N = val.size();
      m.reshape(N, N);
      for(int j = 0; j < N; j++)
      {
        for(int k = 0; k < N; k++)
        {
          m[j][k] = dot(val[j], val[k]);
        }
      }
    };
    m_integrator.integral(i, phi_jk, mass);
  }

  /**
   * \brief 空间的质量矩阵
   */
  void mass_matrix(BSMatrix & mass)
  {
    auto NC = m_mesh->number_of_cells(); 
    auto ldof = m_dof->number_of_local_dofs();
    auto gdof = m_dof->number_of_dofs(); 
    mass.resize(gdof);
    for(int i = 0; i < NC; i++)
    {
      Matrix mat;
      mass_matrix(i, mat);
      auto celldof = m_dof->cell_to_dof(i); 
      for(int j = 0; j < ldof; j++)
      {
        for(int k = 0; k < ldof; k++)
        {
          mass[celldof[j]][celldof[k]] = mat[j][k];
        }
      }
    }
  }

  /**
   * \brief 单元上的质量矩阵
   * \param cidx 单元编号
   * \param mass 返回的质量矩阵
   */
  void curl_matrix(int i, Matrix & curl)
  {
    auto ldof = m_dof->number_of_local_dofs();
    std::function<void(const Node&, Matrix&)> phi_jk = [this, &i]
      (const Node & point, Matrix & m)->void
    {
      std::vector<double> val;
      curl_basis(i, point, val);
      int N = val.size();
      m.reshape(N, N);
      for(int j = 0; j < N; j++)
      {
        for(int k = 0; k < N; k++)
        {
          m[j][k] = val[j]*val[k];
        }
      }
    };
    m_integrator.integral(i, phi_jk, curl);
  }

  /**
   * \brief 空间的质量矩阵
   */
  void curl_matrix(BSMatrix & curl)
  {
    auto NC = m_mesh->number_of_cells(); 
    auto ldof = m_dof->number_of_local_dofs();
    auto gdof = m_dof->number_of_dofs(); 
    curl.resize(gdof);
    for(int i = 0; i < NC; i++)
    {
      Matrix mat;
      curl_matrix(i, mat);
      auto celldof = m_dof->cell_to_dof(i); 
      for(int j = 0; j < ldof; j++)
      {
        for(int k = 0; k < ldof; k++)
        {
          curl[celldof[j]][celldof[k]] = mat[j][k];
        }
      }
    }
  }

  const std::shared_ptr<FKNDof> get_dof()
  {
    return m_dof;
  }

private:
  int m_p;
  Quadrature m_integrator; /**< 积分子 */
  std::shared_ptr<FKNDof> m_dof; /**< 自由度管理 */
  std::shared_ptr<Mesh> m_mesh; /**< 网格 */
  std::shared_ptr<SMSpace> m_smspace; /**< 缩放单项式空间 */
  std::vector<Matrix> m_coeff; /**< 每个单元上基函数的系数 */
};// end of FirstKindNedecFiniteElementSpace2d

}//end of FunctionSpace
}//end of OF

