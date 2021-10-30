
#include <time.h>
#include <initializer_list>
#include <iostream>

#include "TestMacro.h"
#include "common/Permutation.h"

void test_permutation()
{
  OF::Common::Permutation idx({0, 2, 3, 1});
  std::cout<< idx.to_number() <<std::endl;
  std::vector<OF::Common::Permutation> idx1;
  OF::Common::Permutation::gen_all_permutation(4, idx1);

  std::cout<< " idx1 = " <<std::endl;;
  for(int j = 0; j < (int)idx1.size(); j++)
  {
    std::cout<< " [";
    for(int i = 0; i < 4; i++)
    {
      std::cout<< idx1[j][i];
      if(i < 3)
        std::cout<<", ";
    }
    std::cout<< "] : "<< idx1[j].to_number() <<std::endl;
  }
}

int main()
{
  test_permutation();
}
