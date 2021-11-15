
/**
 * \file LagrangeFiniteElementSpace2d.h
 * \brief 二维拉格朗日有限元空间
 * \author 陈春雨
 * \data 2021/11/13
 */
#ifndef LagrangeFiniteElementSpace2d_h
#define LagrangeFiniteElementSpace2d_h

#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <map>

#include <functional>

#include "algebra/linalg.h"
#include "common/Multindex.h"

namespace OF {
namespace FunctionSpace{

template<typename Mesh>
class LFESDof2d
{
public:
  typedef std::array<int, 2> Dof;

  typedef typename OF::Common::MultiIndex MultiIndex;
public:
  LFESDof2d(std::shared_ptr<Mesh> mesh, int p): m_mesh(mesh), m_p(p), m_dpoint(0),
      m_indexs(0)
  {
    get_domain_point();
    MultiIndex::gen_all_index(p, 3, m_indexs);
  }

  /**
   * \brief 获取单元的域点的信息 [维数, 子单形编号, 域点在子单形中的编号]
   */
  void get_domain_point()
  {
    m_dpoint.resize(number_of_local_dofs());
    m_dpoint[0] = {0, 0, 0};

    int num = 1, num0 = 0;
    for(int i = 0; i < m_p-1; i++)
    {
      m_dpoint[num++] = {1, 0, i};
      for(int j = 0; j < i; j++)
      {
        m_dpoint[num++] = {2, 0, num0++}; 
      }
      m_dpoint[num++] = {1, 1, i};
    }
    m_dpoint[num++] = {0, 1, 0};
    for(int i = 0; i < m_p-1; i++)
    {
      m_dpoint[num++] = {1, 2, i};
    }
    m_dpoint[num++] = {0, 2, 0};
  }

  /**
   * \brief 单元的自由度个数
   */
  int number_of_local_dofs()
  {
    return (m_p+1)*(m_p+2)/2;
  }

  /**
   * \brief 每条边内部的自由度个数
   */
  int number_of_edge_local_dofs()
  {
    return m_p>0?m_p-1:0;
  }

  /**
   * \brief 每个单元内部的自由度个数
   */
  int number_of_cell_local_dofs()
  {
    return m_p>1?(m_p-2)*(m_p-1)/2:0;
  }

  /**
   * \brief 网格所有自由度个数
   */
  int number_of_dofs()
  {
    int NC = m_mesh->number_of_cells();
    int NE = m_mesh->number_of_edges();
    int NN = m_mesh->number_of_nodes();

    int NDE = number_of_edge_local_dofs();
    int NDC = number_of_cell_local_dofs();
    return NC*NDC + NE*NDE + NN;
  }

  /**
   * \brief 第 i 条边的自由度的编号
   */
  std::vector<int> edge_to_dof(int i)
  {
    int NDE = number_of_edge_local_dofs();
    int NN = m_mesh->number_of_nodes();
    const auto & edge = m_mesh->edge(i);

    std::vector<int> e2d(m_p+1, 0);
    e2d[0] = edge[0];
    e2d[m_p] = edge[1];
    for(int j = 0; j < NDE; j++)
      e2d[j+1] = NN+i*NDE+j;
    return e2d;
  }

  /**
   * \brief 第 i 个单元的自由度的编号
   */
  std::vector<int> cell_to_dof(int i)
  {
    int NE = m_mesh->number_of_edges();
    int NN = m_mesh->number_of_nodes();
    int NDE = number_of_edge_local_dofs();
    int NDC = number_of_cell_local_dofs();

    int ldof = number_of_local_dofs();
    std::vector<int> c2d(ldof, 0);

    const auto & edge = m_mesh->edges();
    const auto & cell = m_mesh->cell(i);
    const auto & c2e = m_mesh->cell_to_edge(i);

    int idx[3] = {2, 1, 0}; /**< 网格局部边的编号与递增序列编号的索引 TODO */
    bool eflag[3]; /**< 三条边的方向是否正确 */
    eflag[0] = edge[c2e[idx[0]]][0] == cell[0];
    eflag[1] = edge[c2e[idx[1]]][0] == cell[0];
    eflag[2] = edge[c2e[idx[2]]][0] == cell[1];

    for(int j = 0; j < ldof; j++)/**< 循环所有域点 */
    {
      if(m_dpoint[j][0]==0) /**< 0维 域点 */
      {
        c2d[j] = cell[m_dpoint[j][1]];
      }
      else if(m_dpoint[j][0]==1) /**< 1维 域点 */
      {
        int e = c2e[idx[m_dpoint[j][1]]];
        auto e2d = edge_to_dof(e);
        if(eflag[m_dpoint[j][1]])
          c2d[j] = e2d[m_dpoint[j][2]+1];
        else
          c2d[j] = e2d[m_p-m_dpoint[j][2]-1];
      }
      else /**< 2维 域点 */
      {
        c2d[j] = NN + NE*NDE + i*NDC + m_dpoint[j][2];
      }
    }
    return c2d;
  }

private:
  int m_p;
  std::shared_ptr<Mesh> m_mesh;
  std::vector<std::array<int, 3>> m_dpoint;
  std::vector<MultiIndex> m_indexs;/**< 多重指标 */
};// end of FKNDof2d


/**
 * \brief 第一类棱元空间类
 * \param Mesh 空间的网格类
 * \param AK 代数类
 * \param Quadrature 网格上的积分类, 要有单元积分和边积分功能
 */
template<typename Mesh, typename AK, typename Quadrature>
class LagrangeFiniteElementSpace2d
{
public:
  typedef LFESDof2d<Mesh> LFESDof;

  typedef typename Mesh::F F;
  typedef typename Mesh::I I;
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;

  typedef typename Mesh::Vector_2 Vector_2; 

  typedef typename AK::Matrix Matrix;
  typedef typename std::vector<std::map<I, F> > BSMatrix;

  typedef typename OF::Common::MultiIndex MultiIndex;

public:

  /**
   * \brief 有限元空间构造函数
   * \param mesh 空间的网格
   * \param p 空间次数
   * \todo 积分子有两个
   */
  LagrangeFiniteElementSpace2d(std::shared_ptr<Mesh> mesh, int p): 
    m_p(p), m_mesh(mesh), m_integrator(mesh, p+3)
  {
    m_dof = std::make_shared<LFESDof>(mesh, p);
  }

  template<class Container>
  void BMatrix(const Container & bc, std::vector<std::array<double, 3>> & B)
  {
    const int p = m_p;
    /** 计算 B */
    B[0] = {1, 1, 1};
    for(int i = 0; i < p; i++)/**< 赋值 */
    {
      for(int j = 0; j < 3; j++)
        B[i+1][j] = p*bc[j] - i;
    }

    for(int i = 0; i < p; i++)/**< 连乘 */
    {
      for(int j = 0; j < 3; j++)
        B[i+1][j] *= B[i][j];
    }

    int P = 1;
    for(int i = 0; i < p; i++)/**< 乘 P */
    {
      P *= i+1;
      for(int j = 0; j < 3; j++)
        B[i+1][j] /= P;
    }
  }

  template<class Container>
  void nabla_BMatrix(int c, const Container & bc, std::vector<std::array<Vector, 3>> & NB)
  {
    const int p = m_p;
    /** 计算 D^i */
    Matrix D0(p, p);
    Matrix D1(p, p);
    Matrix D2(p, p);
    for(int i = 0; i < p; i++)/**< 赋值 */
    {
      for(int j = 0; j < p; j++)
      {
        D0[i][j] = p*bc[0]-i;
        D1[i][j] = p*bc[1]-i;
        D2[i][j] = p*bc[2]-i;
      }
      D0[i][i] = p;
      D1[i][i] = p;
      D2[i][i] = p;
    }

    for(int i = 1; i < p; i++)/**< 连乘 */
    {
      for(int j = 0; j < p; j++)
      {
        D0[i][j] *= D0[i-1][j];
        D1[i][j] *= D1[i-1][j];
        D2[i][j] *= D2[i-1][j];
      }
    }

    /** 计算 D */
    std::vector<std::array<double, 3> > D(p+1, {0, 0, 0});
    for(int i = 0; i < p; i++)
    {
      for(int j = 0; j < i+1; j++)
      {
        D[i+1][0] += D0[i][j];
        D[i+1][1] += D1[i][j];
        D[i+1][2] += D2[i][j];
      }
    }

    /** 计算 \nabla_lambda */
    Vector_2 v[3];
    Vector GLambda[3];
    v[0] = Vector_2{-1.0, -1.0};
    v[1] = Vector_2{1.0, 0.0};
    v[2] = Vector_2{0.0, 1.0};
    std::vector<std::array<int, 2>> J; 
    m_mesh->jacobi_matrix(c, bc, J);
    if(m_mesh->geo_dimension()==2)
    {
      double det_J = J[0][0]*J[1][1] - J[0][1]*J[1][0];
      for(int i = 0; i < 3; i++)
      {
        GLambda[i][0] = (J[1][1] * v[i][0] - J[1][0] * v[i][1])/det_J;
        GLambda[i][1] = (-J[0][1] * v[i][0] + J[0][0] * v[i][1])/det_J;
      }
    }

    int P = 1;
    for(int i = 0; i < p; i++)/**< 乘 P */
    {
      P *= i+1;
      for(int j = 0; j < 3; j++)
      {
        NB[i+1][j] = D[i+1][j]*GLambda[j]/P;
      }
    }
  }

  /**
   * \brief 计算 Lagrange 基函数函数在重心坐标为 bc 的点处的值.
   * 前三个形函数是重心坐标函数其数值就是 bc[0], bc[1], bc[2].
   */
  template<class Container>
  void basis(int i, const Container & bc, std::vector<double> & val)
  {
    const int p = m_p;
    int ldof = m_dof->number_of_local_dofs();
    /** 计算 B */
    std::vector<std::array<double, 3>> B(p+1, {1, 1, 1});
    BMatrix(bc, B);

    /** 计算基函数值 val */
    MultiIndex midx(p, 3);
    val.resize(ldof);
    for(int i = 0; i < ldof; i++)
    {
      val[i] = B[midx[0]][0]*B[midx[1]][1]*B[midx[2]][2];
      midx.nex_multi_index();
    }
  }

  /**
   * \brief 计算 Lagrange 形函数在重心坐标为 bc 的点处的值.
   * 前三个形函数是重心坐标函数其数值就是 bc[0], bc[1], bc[2].
   */
  template<class Container>
  void grad_basis(int c, const Container & bc, std::vector<Vector> & val)
  {
    const int p = m_p;
    /** 计算 B */
    std::vector<std::array<double, 3>> B(p+1, {1, 1, 1});
    BMatrix(bc, B);

    /** 计算 \nabla B */
    std::vector<std::array<Vector, 3>> NB(p+1);
    nabla_BMatrix(c, bc, NB);

    /** 计算基函数梯度值 val */
    MultiIndex midx(p, 3);
    val.resize((p+1)*(p+2)/2);
    for(int i = 0; i < (int)val.size(); i++)
    {
      val[i] = B[midx[1]][1]*B[midx[2]][2]*NB[midx[0]][0] + 
        B[midx[0]][0]*B[midx[2]][2]*NB[midx[1]][1] + 
        B[midx[0]][0]*B[midx[1]][1]*NB[midx[2]][2];
      midx.nex_multi_index();
    }
  }

  /**
   * \brief 单元上的基函数的质量矩阵
   * \param cidx 单元编号
   * \param mass 返回的质量矩阵
   */
  void mass_matrix(int i, Matrix & mass)
  {
    auto ldof = m_dof->number_of_local_dofs();

    std::function<void(const std::array<double, 3>&, Matrix&)> phi_jk = [this, &i]
      (const std::array<double, 3> & bc, Matrix & m)->void
    {
      std::vector<double> val;
      basis(i, bc, val);
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
    m_integrator.integral(i, phi_jk, mass);
  }

  /**
   * \brief 单元上的基函数的刚度矩阵
   * \param cidx 单元编号
   * \param mass 返回的刚度矩阵
   */
  void stiff_matrix(int i, Matrix & stiff)
  {
    auto ldof = m_dof->number_of_local_dofs();

    std::function<void(const std::array<double, 3>&, Matrix&)> phi_jk = [this, &i]
      (const std::array<double, 3> & bc, Matrix & m)->void
    {
      std::vector<Vector> val;
      grad_basis(i, bc, val);
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
    m_integrator.integral(i, phi_jk, stiff);
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
          mass[celldof[j]][celldof[k]] += mat[j][k];
        }
      }
    }
  }

  /**
   * \brief 空间的刚度矩阵
   */
  void stiff_matrix(BSMatrix & stiff)
  {
    auto NC = m_mesh->number_of_cells(); 
    auto ldof = m_dof->number_of_local_dofs();
    auto gdof = m_dof->number_of_dofs(); 
    stiff.resize(gdof);
    for(int i = 0; i < NC; i++)
    {
      Matrix mat;
      stiff_matrix(i, mat);
      auto celldof = m_dof->cell_to_dof(i); 
      for(int j = 0; j < ldof; j++)
      {
        for(int k = 0; k < ldof; k++)
        {
          stiff[celldof[j]][celldof[k]] += mat[j][k];
        }
      }
    }
  }

  /**
   * \brief 计算右端向量 \f[((f, v_0), (f, v_1), ..., (f, v_{n-1}) \f]
   */
  void source_vector(std::function<double(const Node &)> & f, std::vector<double> & vec)
  {
    auto ldof = m_dof->number_of_local_dofs();
    auto gdof = m_dof->number_of_dofs();
    auto NC = m_mesh->number_of_cells();
    vec.resize(gdof, 0);

    for(int i = 0; i < NC; i++)
    {
      std::function<void(const std::array<double, 3>&, Matrix&)> f_phi_i = [this, &i, f]
        (const std::array<double, 3>  & bc, Matrix & m)->void
      {
        const auto & p0 = m_mesh->node(m_mesh->cell(i)[0]);
        const auto & p1 = m_mesh->node(m_mesh->cell(i)[1]);
        const auto & p2 = m_mesh->node(m_mesh->cell(i)[2]);
        double fval = f(bc[0]*p0 + bc[1]*p1 + bc[2]*p2);

        std::vector<double> pval;
        basis(i, bc, pval);

        int N = pval.size();
        m.reshape(N, 1);
        for(int j = 0; j < N; j++)
        {
          m[j][0] = fval*pval[j];
        }
      };
      Matrix mat;
      m_integrator.integral(i, f_phi_i, mat);
      auto c2d = m_dof->cell_to_dof(i);
      for(int j = 0; j < mat.shape[0]; j++)
      {
        vec[c2d[j]] += mat[j][0];
      }
    }
  }

  const std::shared_ptr<LFESDof> get_dof()
  {
    return m_dof;
  }

private:
  int m_p;
  Quadrature m_integrator; /**< 积分子 */
  std::shared_ptr<LFESDof> m_dof; /**< 自由度管理 */
  std::shared_ptr<Mesh> m_mesh; /**< 网格 */
  std::vector<Matrix> m_coeff; /**< 每个单元上基函数的系数 */
};// end of LagrangeFiniteElementSpace2d

}//end of FunctionSpace
}//end of OF
#endif // end of LagrangeFiniteElementSpace2d_h
