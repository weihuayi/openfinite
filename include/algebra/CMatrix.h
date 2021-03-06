#ifndef CMatrix_h
#define CMatrix_h

#include <iostream>
#include <algorithm>

namespace OF {
namespace AlgebraObject {

template<typename F=double, typename I=int>
struct CMatrix
{
    F * data;
    I shape[2];
    I size;

    CMatrix(I nr, I nc, F val=0.0)
    {
        shape[0] = nr;
        shape[1] = nc;
        size = nr*nc;
        data = new F[size];
        std::fill_n(this->data, size, val);
    }

    ~CMatrix()
    {
        delete [] data;
    }

    void fill(const F val)
    {
        std::fill_n(data, size, val);
    }

    F * operator[](const I i) 
    {
        return &data[i];
    }

    F & operator()(const I i, const I j)
    {
        return data[i*shape[1] + j];
    }

    const F & operator()(const I i, const I j) const
    {
        return data[i*shape[1] + j];
    }

};

template<typename F, typename I>
std::ostream& operator << (std::ostream & os, const CMatrix<F, I> & m)
{
    std::cout << "Matrix("<< m.shape[0] << ","
        << m.shape[1] << "):" << std::endl;
    for(I i = 0; i < m.shape[0]; i ++)
    {
        for(I j = 0; j < m.shape[1]; j++)
        {
            os << m(i, j) << " "; 
        }
        os << "\n";
    }
    os << "\n";
    return os;
}

} // end of namespace AlgebraObject

} // end of namespace OF
#endif // end of CMatrix_h
