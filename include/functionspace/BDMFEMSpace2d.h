
/**
 * \file BDMSpace2d.h
 * \brief 二维 BDM 有限元空间
 * \author 陈春雨
 * \data 2021/11/13
 */
#ifndef BDMSpace2d_h
#define BDMSpace2d_h

#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <map>

#include <functional>

#include "algebra/linalg.h"
#include "common/Multindex.h"
#include "functionspace/ScaledMonomialSpace2d.h"
#include "functionspace/FEMFunction.h"

namespace OF {
namespace FunctionSpace{


/**
 * \brief BDM 有限元的自由度管理
 */
template<typename Mesh>
class BDMDof2d
{
public:
  typedef std::array<int, 2> Dof;
  typedef typename OF::Common::MultiIndex MultiIndex;
  typedef typename Mesh::Vector_2 Vector_2;

public:
  BDMDof2d(const std::shared_ptr<Mesh> mesh, int p): m_mesh(mesh), m_p(p), m_dpoint(0),
      m_indexs(0)
  {
    get_domain_point_and_frame();
    MultiIndex::gen_all_index(p, 3, m_indexs);
  }

  /**
   * \brief 1. 获取单元的域点的信息 [维数, 子单形编号, 域点在子单形中的编号]
   *        2. 获取每个自由度在标准单元的方向.
   *  标准单元:
   *  2
   *  |\
   *  | \
   *  |  \
   *  0---1
   */
  void get_domain_point_and_frame()
  {
    int ND = number_of_local_dofs()/2;
    m_dpoint.resize(ND*2);
    m_direction.resize(ND*2);
    m_tan_flag.resize(ND*2, true);

    m_dpoint[0] = {1, 0, 0};
    m_dpoint[ND] = {1, 1, 0};
    m_direction[0] = Vector_2{0.0, -1.0};
    m_direction[ND] = Vector_2{-1.0, 0.0};
    m_tan_flag[0] = false;
    m_tan_flag[ND] = false;

    int num = 1, num0 = 0;
    double sqrt2 = std::sqrt(2)/2;
    for(int i = 0; i < m_p-1; i++)
    {
      m_dpoint[num] = {1, 0, i+1};
      m_dpoint[num+ND] = {2, 0, num0++};
      m_direction[num] = Vector_2{0.0, -1.0};
      m_direction[num+ND] = Vector_2{1.0, 0.0};
      m_tan_flag[num] = false;
      num++;
      for(int j = 0; j < i; j++)
      {
        m_dpoint[num] = {2, 0, num0++}; 
        m_dpoint[num+ND] = {2, 0, num0++}; 
        m_direction[num] = Vector_2{1.0, 0.0};
        m_direction[num+ND] = Vector_2{0.0, 1.0};
        num++;
      }
      m_dpoint[num] = {1, 1, i+1};
      m_dpoint[num+ND] = {2, 0, num0++};
      m_direction[num] = Vector_2{-1.0, 0.0};
      m_direction[num+ND] = Vector_2{0.0, -1.0};
      m_tan_flag[num] = false;
      num++;
    }

    m_dpoint[num] = {1, 0, m_p};
    m_dpoint[num+ND] = {1, 2, 0};
    m_direction[num] = Vector_2{0.0, -1.0};
    m_direction[num+ND] = Vector_2{sqrt2, sqrt2};
    m_tan_flag[num] = false;
    m_tan_flag[num+ND] = false;
    num++;
    for(int i = 1; i < m_p; i++)
    {
      m_dpoint[num] = {1, 2, i};
      m_dpoint[num+ND] = {2, 0, num0++};
      m_direction[num] = Vector_2{sqrt2, sqrt2};
      m_direction[num+ND] = Vector_2{-sqrt2, sqrt2};
      m_tan_flag[num] = false;
      num++;
    }
    m_dpoint[num] = {1, 1, m_p};
    m_dpoint[num+ND] = {1, 2, m_p};
    m_direction[num] = Vector_2{-1.0, 0.0};
    m_direction[num+ND] = Vector_2{sqrt2, sqrt2};
    m_tan_flag[num] = false;
    m_tan_flag[num+ND] = false;
  }

  /**
   * \brief 单元的自由度个数
   */
  int number_of_local_dofs()
  {
    return (m_p+1)*(m_p+2);
  }

  /**
   * \brief 每条边自由度个数
   */
  int number_of_edge_local_dofs()
  {
    return m_p + 1;
  }

  /**
   * \brief 每个单元内部的自由度个数
   */
  int number_of_cell_local_dofs()
  {
    return number_of_local_dofs() - number_of_edge_local_dofs()*3;
  }

  /**
   * \brief 网格所有自由度个数
   */
  int number_of_dofs()
  {
    int NC = m_mesh->number_of_cells();
    int NE = m_mesh->number_of_edges();

    int NDE = number_of_edge_local_dofs();
    int NDC = number_of_cell_local_dofs();
    return NC*NDC + NE*NDE;
  }

  /**
   * \brief 第 i 条边的自由度的编号
   */
  std::vector<int> edge_to_dof(int i)
  {
    int NDE = number_of_edge_local_dofs();
    std::vector<int> e2d(m_p+1, 0);
    for(int j = 0; j < NDE; j++)
      e2d[j] = i*NDE+j;
    return e2d;
  }

  /**
   * \brief 第 i 个单元的自由度的编号
   */
  std::vector<int> cell_to_dof(int i)
  {
    int NE = m_mesh->number_of_edges();
    int NDE = number_of_edge_local_dofs();
    int NDC = number_of_cell_local_dofs();

    int ldof = number_of_local_dofs();
    std::vector<int> c2d(ldof, 0);

    std::array<int, 3> c2v;
    m_mesh->cell_to_vertex(i, c2v);
    const auto & c2e = m_mesh->cell_to_edge(i);

    std::vector<int> edge2dof[3];
    for(int i = 0; i < 3; i++)
      edge2dof[i] = edge_to_dof(c2e[i]);

    auto & idx = m_mesh->m_localedgeindex; /**< 网格局部边的编号与递增序列编号的索引 TODO */

    bool eflag[3]; /**< 三条边的方向是否正确 */
    const auto & edge = m_mesh->edges();
    eflag[0] = edge[c2e[idx[0]]][0] == c2v[0];
    eflag[1] = edge[c2e[idx[1]]][0] == c2v[0];
    eflag[2] = edge[c2e[idx[2]]][0] == c2v[1];

    for(int j = 0; j < ldof; j++)/**< 循环所有自由度 */
    {
      if(m_dpoint[j][0]==1) /**< 1维 域点 */
      {
        const auto & e2d = edge2dof[idx[m_dpoint[j][1]]];
        if(eflag[m_dpoint[j][1]])
          c2d[j] = e2d[m_dpoint[j][2]];
        else
          c2d[j] = e2d[m_p-m_dpoint[j][2]];
      }
      else /**< 2维 域点 */
      {
        c2d[j] = NE*NDE + i*NDC + m_dpoint[j][2];
      }
    }
    return c2d;
  }

  void is_boundary_dof(std::vector<bool> & isbdof)
  {
    int gdof = number_of_dofs();
    isbdof.resize(gdof, false);
    int NE = m_mesh->number_of_edges();
    for(int i = 0; i < NE; i++)
    {
      if(m_mesh->is_boundary_edge(i))
      {
        auto e2d = edge_to_dof(i); 
        for(auto d : e2d)
        {
          isbdof[d] = true;
        }
      }
    }
  }

  std::vector<int> inverse_dof(int i)
  {
    int ldof = number_of_local_dofs();

    std::array<int, 3> c2v;
    m_mesh->cell_to_vertex(i, c2v);
    const auto & c2e = m_mesh->cell_to_edge(i);

    auto & idx = m_mesh->m_localedgeindex; /**< 网格局部边的编号与递增序列编号的索引 TODO */
    bool eflag[3]; /**< 三条边的方向是否正确 */
    const auto & edge = m_mesh->edges();
    eflag[0] = edge[c2e[idx[0]]][0] == c2v[0];
    eflag[1] = edge[c2e[idx[1]]][0] == c2v[2];
    eflag[2] = edge[c2e[idx[2]]][0] == c2v[1];
    std::cout<< eflag[0] << " " << eflag[1] << " " << eflag[2] <<std::endl;

    std::vector<int> inverse;
    for(int j = 0; j < ldof; j++)/**< 循环所有自由度 */
    {
      if(m_dpoint[j][0]==1) /**< 1维 域点 */
      {
        if(!eflag[m_dpoint[j][1]])
          inverse.push_back(j);
      }
    }
    for(int i = 0; i < (int)inverse.size(); i++)
    {
      std::cout<< "iii" << inverse[i] <<std::endl;
    }
    return inverse;
  }

  std::vector<Vector_2> & get_direction()
  {
    return m_direction;
  }

  std::vector<MultiIndex> & get_multindex()
  {
    return m_indexs;
  }

  std::vector<bool> & get_tan_flag()
  {
    return m_tan_flag;
  }

private:
  int m_p;
  std::shared_ptr<Mesh> m_mesh;
  std::vector<std::array<int, 3>> m_dpoint;
  std::vector<MultiIndex> m_indexs;/**< 多重指标 */
  std::vector<Vector_2> m_direction; /** 标准单元中自由度的方向 */
  std::vector<bool> m_tan_flag;
};// end of BDMDof2d


/**
 * \brief 第一类棱元空间类
 * \param Mesh 空间的网格类
 * \param AK 代数类
 * \param Quadrature 网格上的积分类, 要有单元积分和边积分功能
 */
template<typename Mesh, typename AK, typename Quadrature>
class BDMSpace2d
{
public:
  typedef BDMDof2d<Mesh> BDMDof;

  typedef typename Mesh::F F;
  typedef typename Mesh::I I;
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;

  typedef typename Mesh::Vector_2 Vector_2; 

  typedef typename AK::Matrix Matrix;
  typedef typename std::vector<std::map<I, F> > BSMatrix;

  typedef OF::Common::MultiIndex MultiIndex;
  typedef FEMFunction<BDMSpace2d, Vector> Function;

public:

  /**
   * \brief 有限元空间构造函数
   * \param mesh 空间的网格
   * \param p 空间次数
   * \todo 积分子有两个
   */
  BDMSpace2d(std::shared_ptr<Mesh> mesh, int p): 
    m_p(p), m_mesh(mesh), m_integrator(mesh, p+3)
  {
    m_dof = std::make_shared<BDMDof>(mesh, p);
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
    std::vector<Vector_2> v(3);
    std::vector<Vector> GLambda(3);
    v[0] = Vector_2{-1.0, -1.0};
    v[1] = Vector_2{1.0, 0.0};
    v[2] = Vector_2{0.0, 1.0};
    m_mesh->push_dual_vector(c, bc, v, GLambda);

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
  void basis(int c, const Container & bc, std::vector<Vector> & val)
  {
    /** 计算 B */
    std::vector<std::array<double, 3>> B(m_p+1, {1, 1, 1});
    BMatrix(bc, B);

    m_mesh->push_frame(c, bc, m_dof->get_direction(), val, m_dof->get_tan_flag());
    auto invdof = m_dof->inverse_dof(c);
    for(auto didx : invdof)
      val[didx] *=  -1;

    /** 计算基函数值 val */
    int hldof = m_dof->number_of_local_dofs()/2;
    for(int i = 0; i < hldof; i++)
    {
      auto & midx = m_dof->get_multindex()[i];
      val[i] *= B[midx[0]][0]*B[midx[1]][1]*B[midx[2]][2];
      val[hldof+i] *=  B[midx[0]][0]*B[midx[1]][1]*B[midx[2]][2];
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
      std::vector<Vector> val;
      basis(i, bc, val);
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
  template<typename Container>
  void source_vector(std::function<Vector(const Node &)> & f, Container & vec)
  {
    auto ldof = m_dof->number_of_local_dofs();
    auto gdof = m_dof->number_of_dofs();
    auto NC = m_mesh->number_of_cells();
    if(vec.size() != gdof)
      vec.resize(gdof, 0);

    for(int i = 0; i < NC; i++)
    {
      std::function<void(const std::array<double, 3>&, Matrix&)> f_phi_i = [this, &i, f]
        (const std::array<double, 3>  & bc, Matrix & m)->void
      {
        auto po = m_mesh->bc_to_point(i, bc);
        Vector fval = f(po);

        std::vector<Vector> pval;
        basis(i, bc, pval);

        int N = pval.size();
        m.reshape(N, 1);
        for(int j = 0; j < N; j++)
        {
          m[j][0] = dot(fval, pval[j]);
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

  Function function()
  {
    return Function(this);
  }

  double L2_error(Function & Fun, std::function<Vector(const Node &)> & fun)
  {
    double val = 0.0;
    int ldof = m_dof->number_of_local_dofs();
    int NC = m_mesh->number_of_cells();
    for(int i = 0; i < NC; i++)
    {
      std::function<double(const std::array<double, 3>&)> L2_err = [ldof, this, i, Fun, fun]
        (const std::array<double, 3>  & bc)->double
      {
        auto po = m_mesh->bc_to_point(i, bc);
        return (Fun(i, bc)-fun(po)).squared_length();
      };
      val += m_integrator.integral(i, L2_err);
      std::cout<< i << " : " << val <<std::endl;
    }
    return std::sqrt(val);
  }

  Function inteprate(std::function<Vector(const Node &)> & f)
  {
    Function F(this);
    auto & val = F.get_data();

    int gdof = m_dof->number_of_dofs();
    int ldof = m_dof->number_of_local_dofs();
    val.resize(gdof, 0);
    int NC = m_mesh->number_of_cells();

    std::vector<MultiIndex> midxs;
    MultiIndex::gen_all_index(m_p, 3, midxs);
    std::array<double, 3> bc = {0.0, 0.0, 0.0};

    for(int i = 0; i < NC; i++)
    {
      auto c2d = m_dof->cell_to_dof(i);
      std::vector<Vector_2> v0(2);
      std::vector<Vector> v1(2);
      std::vector<bool> tf(2);
      for(int j = 0; j < ldof/2; j++)
      {
        auto & idx = midxs[j];
        bc[0] = (double)idx[0]/m_p;
        bc[1] = (double)idx[1]/m_p;
        bc[2] = (double)idx[2]/m_p;
        v0[0] = m_dof->get_direction()[j];
        v0[1] = m_dof->get_direction()[j+ldof/2];
        tf[0] = m_dof->get_tan_flag()[j];
        tf[1] = m_dof->get_tan_flag()[j+ldof/2];

        m_mesh->push_frame(i, bc, v0, v1, tf);
        auto po = m_mesh->bc_to_point(i, bc);
        auto fval = f(po);

        double det = v1[0][0]*v1[1][1] - v1[0][1]*v1[1][0]; 
        val[c2d[j]] = (v1[1][1]*fval[0] - v1[1][0]*fval[1])/det;
        val[c2d[j+ldof/2]] = (-v1[0][1]*fval[0] + v1[0][0]*fval[1])/det;
      }
    }
    return F;
  }

  const std::shared_ptr<BDMDof> get_dof()
  {
    return m_dof;
  }

private:
  int m_p;
  Quadrature m_integrator; /**< 积分子 */
  std::shared_ptr<BDMDof> m_dof; /**< 自由度管理 */
  std::shared_ptr<Mesh> m_mesh; /**< 网格 */
  std::vector<Matrix> m_coeff; /**< 每个单元上基函数的系数 */
};// end of BDMSpace2d

}//end of FunctionSpace
}//end of OF
#endif // end of BDMSpace2d_h
