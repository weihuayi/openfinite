#include <iostream>
#include <map>

#include "core/CuthillMckee.h"
typedef OF::Core::BSMatrix BSMatrix;

void get_A(BSMatrix & A)
{
 // A.resize(21);
 // for(int i = 0; i < 21; i ++)
 // {
 //   A[i][i] = 1.0;
 //   if(i%7<6)
 //   {
 //     A[i][i+1] = 1.0;
 //     A[i+1][i] = 1.0;
 //   }
 //   if(i/7<2)
 //   {
 //     A[i][i+7] = 1.0;
 //     A[i+7][i] = 1.0;
 //   }
 // }
 A.resize(6);
 A[0][0] = 1.0;
 A[0][5] = 1.0;
 A[1][1] = 1.0;
 A[2][0] = 1.0;
 A[2][2] = 1.0;
 A[1][2] = 1.0;
 A[2][1] = 1.0;
 A[5][0] = 1.0;
 A[5][5] = 1.0;
 A[0][3] = 1.0;
 A[3][0] = 1.0;
 A[4][0] = 1.0;
 A[0][4] = 1.0;
 A[4][4] = 1.0;
 A[3][3] = 1.0;
}

int main()
{
  BSMatrix A, B; 
  get_A(A);

  std::vector<int> perm;
  A.cuthill_mckee(perm);

  int N = perm.size();
  std::vector<int> rperm(N); /**< perm 的逆排列 */
  for(int i = 0; i < N; i++)
    rperm[perm[i]] = i;

  A.rearrange(rperm, B);

  for(int i = 0; i < N; i++)
  {
    std::cout<< perm[i]<< " ";
  }
  std::cout<< "" <<std::endl;

  std::vector<double> A0, B0; 
  A.todense(A0);
  B.todense(B0);

  std::cout << "A = " << std::endl;
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
    {
      std::cout << A0[i*N+j] << " ";
    }
    std::cout << "" << std::endl;
  }
  std::cout << "B = " << std::endl;

  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
    {
      std::cout << B0[i*N+j] << " ";
    }
    std::cout << "" << std::endl;
  }
}








