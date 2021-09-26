#ifndef ScaledMonomialSpace2d_h
#define ScaledMonomialSpace2d_h

#include <math.h>
#include <memory>
#include <vector>
#include <array>

#include <functional>

namespace OF {
namespace FunctionSpace{

template<typename Mesh>
class SMDof2d
{
public:
  typedef std::array<int, 2> Dof;

public:
  SMDof2d(std::shared_ptr<Mesh> mesh, int p): m_mesh(mesh), m_p(p) {}

  int number_of_local_dofs(int p = -1)
  {
    if(p < 0)
      p = m_p;
    return (p+1)*(p+2)/2;
  }

  int number_of_dofs()
  {
    auto NC = m_mesh->number_of_cells();
    return NC*number_of_local_dofs();
  }

  void cell_to_dof(int idx, std::vector<int> & cell2dof)
  {
    int ldof = number_of_local_dofs();
    cell2dof.resize(ldof);

    int N = idx*ldof;
    for(int i = 0; i < ldof; i++)
    {
      cell2dof[i] = N+i;
    }
  }

  /**
   * \brief 第 i 个单元的第 j 个自由度的编号
   */
  int cell_to_dof(int i, int j)
  {
    return i*number_of_local_dofs()+j;
  }


private:
  int m_p;
  std::shared_ptr<Mesh> m_mesh;
};// end of SMDof2d

template<typename Mesh, typename AK, typename Quadrature>
class ScaledMonomialSpace2d
{
public:
  typedef SMDof2d<Mesh> SMDof;

  typedef typename Mesh::F F;
  typedef typename Mesh::I I;
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;

  typedef typename AK::Matrix Matrix;

public:
  ScaledMonomialSpace2d(std::shared_ptr<Mesh> mesh, int p): 
    m_mesh(mesh), m_p(p), m_integrator(mesh, p+3)
  {
    m_dof = std::make_shared<SMDof>(mesh, p); /**< \todo 为什么放在初始化里不行 */
    mesh->cell_barycenter(m_cellbarycenter);
    mesh->cell_size(m_cellsize);
    mesh->cell_measure(m_cellmeasure);
  }

  void basis(int cidx, const Node & point, std::vector<double> & phi, int p = -1)
  {
    if(p < 0)
      p = m_p;

    int ldof = m_dof->number_of_local_dofs(p);
    phi.resize(ldof);

    const auto & cellbar = m_cellbarycenter[cidx];
    const auto & h = m_cellsize[cidx];
    auto xbar = (point - cellbar)/h;

    phi[0] = 1; // 0 次基函数
    if(p>0)
    {
      phi[1] = xbar[0];
      phi[2] = xbar[1];// 1 次基函数
      int start = 3; // 第 0 个 j 次基函数的编号
      for(int j = 2; j < p+1; j++)
      {
        for(int k = start; k < start+j; k++)
        {
          phi[k] = xbar[0]*phi[k-j]; 
        }
        phi[start+j] = xbar[1]*phi[start-1];
        start += j+1;
      }
    }
  }

  void grad_basis(int cidx, const Node & point, std::vector<Vector> & nphi, int p=-1)
  {
    if(p < 0)
      p = m_p;

    int ldof = m_dof->number_of_local_dofs(p);
    nphi.resize(ldof);

    const auto & cellbar = m_cellbarycenter[cidx];
    const auto & h = m_cellsize[cidx];
    auto xbar = (point - cellbar)/h;

    std::vector<double> phi; //基函数 
    basis(cidx, point, phi);

    nphi[0][0] = 0.0;
    nphi[0][1] = 0.0; // 0 次基函数的偏导
    if(p>0)
    {
      int start = 1; // 第 0 个 j 次基函数的编号
      for(int j = 1; j < p+1; j++)
      {
        nphi[start][1] = 0.0; // phi[-1]
        for(int k = start; k < start+j; k++)
        {
          int c = k - start;
          nphi[k][0] = (j-c)*phi[k-j];
          nphi[k+1][1] = (c+1)*phi[k-j];
        }
        nphi[start+j][0] = 0.0;
        start += j+1; //更新 start
      }
    }

    for(int j = 0; j < ldof; j++)
      nphi[j] /= h;
  }

  void laplace_basis(int cidx, const Node & point, std::vector<double> & lphi)
  {
    int ldof = m_dof->number_of_local_dofs();
    lphi.resize(ldof);

    const auto & cellbar = m_cellbarycenter[cidx];
    const auto & h = m_cellsize[cidx];
    const auto & area = m_cellmeasure[cidx];
    auto xbar = (point - cellbar)/h;

    std::vector<double> phi; //基函数 
    basis(cidx, point, phi);

    lphi[0] = 0.0; // 0 次基函数的偏导
    if(m_p>0)
    {
      lphi[1] = 0.0;
      lphi[2] = 0.0;
      phi[1] = xbar[0];
      phi[2] = xbar[1];
      int start = 3; // 第 0 个 j 次基函数的编号
      for(int j = 2; j < m_p+1; j++)
      {
        for(int k = start; k < start+j-1; k++)
        {
          int c = k - start;
          lphi[k] += (j-c)*(j-c-1)*phi[k-j-j+1];
          lphi[k+2] += (c+2)*(c+1)*phi[k-j-j+1];
        }
        start += j+1; //更新 start
      }
    }

    for(int j = 0; j < ldof; j++)
      lphi[j] /= area;
  }

  /**
   * \brief 边上的基函数: 1, x, x**2, x**3, ..., x**{n-1}
   * \param eidx 边的编号
   * \param point 节点笛卡尔坐标
   * \param phi 返回值
   */
  void edge_basis(int eidx, const Node & point, std::vector<double> & phi)
  {
    int ldof = m_p+1;
    phi.resize(ldof);

    const auto & edgebar = m_mesh->edge_barycenter(eidx);
    const auto & h = m_mesh->edge_measure(eidx);
    auto xbar = dot(point - edgebar, m_mesh->edge_tangent(eidx))/h;
    phi[0] = 1;
    for(int i = 1; i < ldof; i++)
    {
      phi[i] = xbar*phi[i-1];
    }
  }
  
  /**
   * \brief 空间的质量矩阵
   * \todo 将 mass 改为稀疏矩阵
   */
  void mass_matrix(int i, Matrix & mass, int p = -1)
  {
    if(p < 0)
      p = m_p;

    auto ldof = m_dof->number_of_local_dofs(p);
    std::function<void(const Node&, Matrix&)> phi_jk = [this, &i, &p]
      (const Node & point, Matrix & m)->void
    {
      std::vector<F> val;
      basis(i, point, val, p);
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
   * \brief 空间的质量矩阵
   * \todo 将 mass 改为稀疏矩阵
   */
  void mass_matrix(Matrix & mass)
  {
    auto NC = m_mesh->number_of_cells(); 
    auto ldof = m_dof->number_of_local_dofs();
    auto gdof = NC*ldof; 
    mass.reshape(gdof, gdof); //TODO mass 应该是一个稀疏矩阵
    for(int i = 0; i < NC; i++)
    {
      Matrix mat;
      mass_matrix(i, mat);
      for(int j = 0; j < ldof; j++)
      {
        auto jdx = m_dof->cell_to_dof(i, j);
        for(int k = 0; k < ldof; k++)
        {
          auto kdx = m_dof->cell_to_dof(i, k);
          mass[jdx][kdx] = mat[j][k];
        }
      }
    }
  }

  /**
   * \brief 空间的刚度矩阵
   * \todo 将 stiff 改为稀疏矩阵
   */
  void stiff_matrix(Matrix & stiff)
  {
    auto NC = m_mesh->number_of_cells(); 
    auto ldof = m_dof->number_of_local_dofs();
    auto gdof = NC*ldof; 
    stiff.reshape(gdof, gdof); //TODO stiff 应该是一个稀疏矩阵
    for(int i = 0; i < NC; i++)
    {
      std::function<void(const Node&, Matrix&)> phi_jk = [this, i]
        (const Node & point, Matrix & m)->void
      {
        std::vector<Vector> val;
        this->grad_basis(i, point, val);
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

      Matrix mat;
      m_integrator.integral(i, phi_jk, mat);
      for(int j = 0; j < ldof; j++)
      {
        auto jdx = m_dof->cell_to_dof(i, j);
        for(int k = 0; k < ldof; k++)
        {
          auto kdx = m_dof->cell_to_dof(i, k);
          stiff[jdx][kdx] = mat[j][k];
        }
      }
    }
  }

private:
  int m_p;
  Quadrature m_integrator;
  std::vector<Node> m_cellbarycenter;
  std::vector<F> m_cellsize;
  std::vector<F> m_cellmeasure;
  std::shared_ptr<Mesh> m_mesh;
  std::shared_ptr<SMDof> m_dof;
};// end of ScaledMonomialSpace2d

}//end of FunctionSpace
}//end of OF
#endif // end of ScaledMonomialSpace2d_h
