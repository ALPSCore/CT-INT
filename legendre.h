//
// Created by H. Shinaoka on 2015/10/21.
//

#ifndef IMPSOLVER_LEGENDRE_H
#define IMPSOLVER_LEGENDRE_H

#include<complex>
#include<cmath>
#include<vector>
#include<assert.h>

#include "boost/math/special_functions/bessel.hpp"
#include "boost/multi_array.hpp"

#include <alps/numeric/matrix.hpp>


class LegendreTransformer {
  public:
    LegendreTransformer(int n_matsubara, int n_legendre);

  private:
    const int n_matsubara_, n_legendre_;

public:
    const alps::numeric::matrix<std::complex<double> > & Tnl() const;
    void compute_legendre(double x, std::vector<double>& val) const;
    void compute_legendre(const std::vector<double>& xval, boost::multi_array<double,2>& val) const;

private:
    alps::numeric::matrix<std::complex<double> > Tnl_;
    std::vector<double> inv_l_;
};


#endif //IMPSOLVER_LEGENDRE_H
