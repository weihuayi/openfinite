
/**
 * \file CuthillMckee.h
 * \brief 三维拉格朗日有限元空间
 * \author Chunyu Chen
 * \data 2022/2/13
 */
#ifndef CuthillMckee_h
#define CuthillMckee_h

#include <vector>
#include <queue>
#include <algorithm>
#include <map>
#include <iostream>

namespace OF {
namespace Core{

class BSMatrix: public std::vector<std::map<int, double> >
{
public:
  BSMatrix(){};

  void cuthill_mckee(std::vector<int> & perm)
  {
    int N = this->size();
    perm.resize(N, -1);
    std::vector<int> idx(N, 0);
    for(int i = 0; i < N; i++)
      idx[i] = i;

    auto f0 = [this](int a, int b) 
    {
      return this->at(a).size() < this->at(b).size();
    };
    int v = *std::min_element(idx.begin(), idx.end(), f0);
    perm[v] = 0;

    std::vector<int> S;
    std::queue<int> R;
    R.push(v);

    int k = 1;
    while(!R.empty())
    {
      int x = R.front();
      std::cout << "x = " << x << std::endl;
      for(const auto & m : this->at(x))
      {
        if(perm[m.first] < 0)
          S.push_back(m.first);
      }
      std::sort(S.begin(), S.end(), f0);
      for(const auto & y : S)
      {
        std::cout << "y = " << y << std::endl;
        perm[y] = k++;
        R.push(y);
      }
      R.pop();
      S.resize(0);
    }
  }

  void rearrange(const std::vector<int> & perm, BSMatrix & out)
  {
    int N = this->size();
    std::vector<int> rperm(N); /**< perm 的逆排列 */
    for(int i = 0; i < N; i++)
      rperm[perm[i]] = i;

    out.resize(N);
    for(int i = 0; i < N; i++)
    {
      for(auto val : this->at(perm[i]))
      {
        out[i][rperm[val.first]] = val.second;
      }
    }
  }

  void todense(std::vector<double> & v)
  {
    int N = this->size();
    v.resize(N*N);
    for(int i = 0; i < N; i++)
    {
      for(auto val : this->at(i))
      {
        int idx = val.first + i*N;
        v[idx] = val.second;
      }
    }
  }
};

}
}




#endif // end of CuthillMckee_h
