//
// Created by Hiroshi Shinaoka on 2015/10/15.
//

#ifndef IMPSOLVER_UTIL_H
#define IMPSOLVER_UTIL_H

#include <vector>
#include <complex>

#include <alps/numeric/matrix.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>

template<typename T> T mycast(std::complex<double> val);
template<typename T> T myconj(T val);

//typedef alps::numeric::matrix<double> dense_matrix;
//typedef alps::numeric::matrix<std::complex<double> > complex_dense_matrix;
//
//template<class MAT, class SCALAR>
//inline void invert(MAT & A, SCALAR & det);

template<class R>
std::vector<size_t> pickup_a_few_numbers(size_t N, size_t n, R& random01) {
    std::vector<size_t> flag(N,0), list(n);

    for (size_t i=0; i<n; ++i) {
        size_t itmp = 0;
        while(true) {
            itmp = std::min((size_t)(random01()*N),N-1);
            if (flag[itmp]==0) {
                break;
            }
        }
        list[i] = itmp;
        flag[itmp] = 1;
    }
    return list;
}

namespace alps {
    namespace numeric {

        template<class Matrix>
        typename Matrix::value_type determinant(Matrix M) {
            std::vector<int> ipiv(num_rows(M));
            const size_t N = num_rows(M);

            if (N==0) {
                return 1.0;
            }

            int info = boost::numeric::bindings::lapack::getrf(M, ipiv);
            if (info != 0)
                throw std::runtime_error("Error in GETRF !");

            double det = 1.0;
            for (size_t i=0; i<N; ++i) {
                det *= M(i,i);
            }
            int p = 0;
            for (size_t i=0; i<N-1; ++i) {
                if (ipiv[i] != i+1) {
                    ++p;
                }
            }
            return (p%2==0 ? det : -det);
        }
    }
}

#endif //IMPSOLVER_UTIL_H
