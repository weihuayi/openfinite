#ifndef FEMFunction_h
#define FEMFunction_h

#include <vector>
#include <array>
namespace OF{
namespace FunctionSpace{

template<typename Space, typename Val>
class FEMFunction
{
public:
  FEMFunction(Space * space):m_space(space), m_data(0){}

  void set_data(std::vector<double> & data)
  {
    m_data = data;
  }

  std::vector<double> & get_data()
  {
    return m_data;
  }

  double & operator[](const int & i)
  {
    return m_data[i];
  }

  Val operator()(int c, const std::array<double, 3> & bc)const
  {
    std::vector<Val> val; 
    m_space->basis(c, bc, val);
    auto dof = m_space->get_dof()->cell_to_dof(c);
    val[0] *= m_data[dof[0]];
    for(int i = 1; i < (int)dof.size(); i++)
    {
      val[i] *= m_data[dof[i]];
      val[0] += val[i];
    }
    return val[0];
  }

  Val value(int i, const std::array<double, 3> & bc)const
  {
    std::vector<Val> val; 
    m_space->basis(i, bc, val);
    auto dof = m_space->get_dof()->cell_to_dof(i);
    val[0] *= m_data[dof[0]];
    for(int i = 1; i < (int)dof.size(); i++)
    {
      val[0] += val[i]*m_data[dof[i]];
    }
    return val[0];
  }

private:
  Space * m_space;
  std::vector<double> m_data;
}; // end of FEMFunction

}//end of FunctionSpace
}//end of OF
#endif // end of FEMFunction_h
