//
// Created by H. Shinaoka on 2015/10/15.
//

#include <vector>
#include <assert.h>
#include "util.h"

template<>
double mycast(std::complex<double> val) {
    return val.real();
}

template<>
std::complex<double> mycast(std::complex<double> val) {
    return val;
}

template<>
double myconj(double val) {
    return val;
}

template<>
std::complex<double> myconj(std::complex<double> val) {
    return std::conj(val);
}

/*
template<>
inline void invert<dense_matrix,double>(dense_matrix & A, double & det) {
    dense_matrix B(A.size1(), A.size1());

    B = boost::numeric::ublas::identity_matrix<double>(A.size1());

    if (A.size1()>0) {
        boost::numeric::ublas::vector<fortran_int_t> ipivot(A.size1());
        boost::numeric::bindings::lapack::gesv(A, ipivot, B);
        swap(A,B);
        det = 1;
        for (int i=0; i<(int)A.size1(); i++) {
            det *= B(i,i);
        }
        det = std::fabs(det);
    }
    else {
        det = 1;
    }
}

template<>
inline void invert<complex_dense_matrix,std::complex<double> >(complex_dense_matrix & A, std::complex<double> & det) {
    complex_dense_matrix B(A.size1(), A.size1());

    B = boost::numeric::ublas::identity_matrix<double>(A.size1());

    if (A.size1()>0) {
        boost::numeric::ublas::vector<fortran_int_t> ipivot(A.size1());
        boost::numeric::bindings::lapack::gesv(A, ipivot, B);
        swap(A,B);
        det = 1.;
        for (int i=0; i<(int)A.size1(); i++) {
            det *= B(i,i);
        }
    }
    else {
        det = 1.;
    }
}
*/
