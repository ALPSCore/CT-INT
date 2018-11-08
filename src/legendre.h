//
// Created by H. Shinaoka on 2015/10/21.
//

#pragma once

#include<complex>
#include<cmath>
#include<vector>
#include<assert.h>

#include "boost/math/special_functions/bessel.hpp"
#include "boost/multi_array.hpp"

#include <Eigen/Core>

namespace alps {
    namespace ctint {

        class LegendreTransformer {
        public:
            LegendreTransformer(int n_matsubara, int n_legendre)
              : n_matsubara_(n_matsubara), n_legendre_(n_legendre), Tnl_(n_matsubara, n_legendre), inv_l_(n_legendre) {

              double sign_tmp = 1.0;
              for (int im = 0; im < n_matsubara_; ++im) {
                std::complex<double> ztmp(0.0, 1.0);
                for (int il = 0; il < n_legendre_; ++il) {
                  Tnl_(im, il) =
                    sign_tmp * ztmp * std::sqrt(2 * il + 1.0) * boost::math::sph_bessel(il, 0.5 * (2 * im + 1) * M_PI);
                  ztmp *= std::complex<double>(0.0, 1.0);
                }
                sign_tmp *= -1;
              }
              sqrt_2l_1.resize(n_legendre);
              sqrt_2l_1[0] = 1.0;
              for (int l = 1; l < n_legendre_; l++) {
                inv_l_[l] = 1.0 / l;
                sqrt_2l_1[l] = std::sqrt(2.0 * l + 1.0);
              }
            };

        private:
            const int n_matsubara_, n_legendre_;

        public:
            const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> &Tnl() const {
              return Tnl_;
            }

            int Nl() const {
              return n_legendre_;
            }

            void
            compute_legendre(double x, std::vector<double> &val) const {
              assert(val.size() >= n_legendre_);
              for (int l = 0; l < n_legendre_; l++) {
                if (l == 0) {
                  val[l] = 1;
                } else if (l == 1) {
                  val[l] = x;
                } else {
                  //val[l] = ((2 * l - 1) * x * val[l-1] - (l - 1) * val[l-2]) / static_cast<double>(l);//l
                  val[l] = ((2 * l - 1) * x * val[l - 1] - (l - 1) * val[l - 2]) * inv_l_[l];//l
                }
              }
            }

            template <typename MULTI_ARRAY>
            void
            compute_legendre(const std::vector<double> &xval,
                                                  MULTI_ARRAY &val) const {
              assert(val.shape()[0] >= n_legendre_);
              const int nx = xval.size();
              for (int l = 0; l < n_legendre_; l++) {
                if (l == 0) {
#pragma clang loop vectorize(enable)
                  for (int ix = 0; ix < nx; ++ix) {
                    val[l][ix] = 1;
                  }
                } else if (l == 1) {
#pragma clang loop vectorize(enable)
                  for (int ix = 0; ix < nx; ++ix) {
                    val[l][ix] = xval[ix];
                  }
                } else {
                  //for (int ix=0; ix<nx; ++ix) {
                  //val[ix][l] = ((2 * l - 1) * xval[ix]*val[ix][l - 1] - (l - 1) * val[ix][l - 2]) * inv_l_[l];//l
                  //}
                  const double inv_l_tmp = inv_l_[l];
#pragma clang loop vectorize(enable)
                  for (int ix = 0; ix < nx; ++ix) {
                    val[l][ix] = ((2 * l - 1) * xval[ix] * val[l - 1][ix] - (l - 1) * val[l - 2][ix]) * inv_l_tmp;//l
                  }
                }
              }
            }

            const std::vector<double>& get_sqrt_2l_1() const {return sqrt_2l_1;}

        private:
            Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> Tnl_;
            std::vector<double> inv_l_, sqrt_2l_1;
        };

    }
}

