//
// Created by H. Shinaoka on 2015/10/21.
//

#include "legendre.h"


LegendreTransformer::LegendreTransformer(int n_matsubara, int n_legendre)
        : n_matsubara_(n_matsubara), n_legendre_(n_legendre), Tnl_(n_matsubara,n_legendre), inv_l_(n_legendre) {

    double sign_tmp = 1.0;
    for (int im=0; im<n_matsubara_; ++im) {
        std::complex<double> ztmp(0.0,1.0);
        for (int il=0; il<n_legendre_; ++ il) {
            Tnl_(im,il) = sign_tmp*ztmp*std::sqrt(2*il+1.0)*boost::math::sph_bessel(il,0.5*(2*im+1)*M_PI);
            ztmp *= std::complex<double>(0.0,1.0);
        }
        sign_tmp *= -1;
    }
    for(int l=1; l<n_legendre_; l++) {
        inv_l_[l] = 1.0/l;
    }
};

void
LegendreTransformer::compute_legendre(double x, std::vector<double> &val) const {
    assert(val.size()>=n_legendre_);
    for(int l=0; l<n_legendre_; l++) {
        if (l == 0) {
            val[l] = 1;
        } else if (l == 1) {
            val[l] = x;
        } else {
            //val[l] = ((2 * l - 1) * x * val[l-1] - (l - 1) * val[l-2]) / static_cast<double>(l);//l
            val[l] = ((2 * l - 1) * x * val[l-1] - (l - 1) * val[l-2])*inv_l_[l];//l
        }
    }
}

const alps::numeric::matrix<std::complex<double> > &LegendreTransformer::Tnl() const {
    return Tnl_;
}

void
LegendreTransformer::compute_legendre(const std::vector<double>& xval, boost::multi_array<double,2> &val) const {
    assert(val.shape()[0]>=n_legendre_);
    const int nx = xval.size();
    for(int l=0; l<n_legendre_; l++) {
        if (l == 0) {
#pragma clang loop vectorize(enable)
            for (int ix=0; ix<nx; ++ix) {
                val[l][ix] = 1;
            }
        } else if (l == 1) {
#pragma clang loop vectorize(enable)
            for (int ix=0; ix<nx; ++ix) {
                val[l][ix] = xval[ix];
            }
        } else {
            //for (int ix=0; ix<nx; ++ix) {
                //val[ix][l] = ((2 * l - 1) * xval[ix]*val[ix][l - 1] - (l - 1) * val[ix][l - 2]) * inv_l_[l];//l
            //}
            const double inv_l_tmp = inv_l_[l];
#pragma clang loop vectorize(enable)
            for (int ix=0; ix<nx; ++ix) {
                val[l][ix] = ((2 * l - 1) * xval[ix]*val[l - 1][ix] - (l - 1) * val[l-2][ix]) * inv_l_tmp;//l
            }
        }
    }
}

