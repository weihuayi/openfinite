
#include <time.h>
#include <initializer_list>
#include <iostream>

#include "TestMacro.h"
#include "common/IncreaseSequence.h"

void test_permutation(int l = 1, int n = 6)
{
  OF::Common::IncreaseSequence idx({0, 2, 3});
  OF::Common::IncreaseSequence idx0({0, 1, 2});
  idx.nex_increaseing_sequence(4, idx0);
  for(int i = 0; i < l; i++)
    std::cout<< idx0[i] <<std::endl;
  std::vector<OF::Common::IncreaseSequence> idx1;
  OF::Common::IncreaseSequence::gen_all_sequence(n, l, idx1);
  std::cout<< " idx1 = " <<std::endl;;
  for(int j = 0; j < (int)idx1.size(); j++)
  {
    std::cout<< " [";
    for(int i = 0; i < l; i++)
    {
      std::cout<< idx1[j][i];
      if(i < l-1)
        std::cout<<", ";
    }
    std::cout<< "] : "<< idx1[j].to_number(n) <<std::endl;
  }
  /*
  std::cout<< "idx0 = " <<std::endl;
  for(int i = 0; i < l; i++)
    std::cout<< idx1[2][i] <<std::endl;
  std::cout<< "idx0num = " <<std::endl;
  std::cout<< idx1[2].to_number(n) <<std::endl;
  */
}

int main()
{
  test_permutation();
}
