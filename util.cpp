//
// Created by H. Shinaoka on 2015/10/15.
//

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
