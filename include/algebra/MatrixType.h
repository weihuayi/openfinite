#ifndef MatrixType_h
#define MatrixType_h

namespace OF {
namespace AlgebraObject {

enum MatrixType  // 读取文件行的类型
{
  F,
  FR,
  FC,
  CSR,
  CSC,
  COO,
  BSR,
  BSC
};

} // end of namespace Algebra

} // end of namespace OF
#endif // end of MatrixType_h
