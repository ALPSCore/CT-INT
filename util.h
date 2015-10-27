//
// Created by Hiroshi Shinaoka on 2015/10/15.
//

#ifndef IMPSOLVER_UTIL_H
#define IMPSOLVER_UTIL_H

#include <vector>
#include <complex>

#include <alps/numeric/matrix.hpp>

template<typename T> T mycast(std::complex<double> val);
template<typename T> T myconj(T val);

//typedef alps::numeric::matrix<double> dense_matrix;
//typedef alps::numeric::matrix<std::complex<double> > complex_dense_matrix;
//
//template<class MAT, class SCALAR>
//inline void invert(MAT & A, SCALAR & det);

#endif //IMPSOLVER_UTIL_H
