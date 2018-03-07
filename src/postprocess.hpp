//
// Created by H. Shinaoka on 2015/12/23.
//

#pragma once

#include <complex>

//#include "interaction_expansion.hpp"
#include "legendre.h"

namespace alps {
    namespace ctint {

        template<class SOLVER_TYPE>
        void evaluate_selfenergy_measurement_legendre(const typename alps::accumulators::result_set &results,
                                                      //typename SOLVER_TYPE::matsubara_green_function_t &green_matsubara_measured,
                                                      //typename SOLVER_TYPE::matsubara_green_function_t &sigma_green_matsubara_measured,
                                                      //typename SOLVER_TYPE::itime_green_function_t &green_itime_measured,
                                                      //const typename SOLVER_TYPE::matsubara_green_function_t &bare_green_matsubara,
                                                      const double &beta,
                                                      int n_site,
                                                      int n_flavors,
                                                      int n_matsubara,
                                                      int n_legendre,
                                                      boost::multi_array<std::complex<double>,4>& Sl,
                                                      boost::multi_array<std::complex<double>,4>& Sw
        ) {
          std::cout << "evaluating self energy measurement: lengendre, real space" << std::endl;
          double sign = results["Sign"].template mean<double>();

          //Legendre expansion utils
          LegendreTransformer legendre_transformer(n_matsubara, n_legendre);
          auto &Tnl = legendre_transformer.Tnl();

          Sl = boost::multi_array<std::complex<double>, 4>(boost::extents[n_legendre][n_site][n_site][n_flavors]);
          Sw = boost::multi_array<std::complex<double>, 4>(boost::extents[n_matsubara][n_site][n_site][n_flavors]);

          //load S_l
          for (int flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
            for (int site1 = 0; site1 < n_site; ++site1) {
              for (int site2 = 0; site2 < n_site; ++site2) {
                std::stringstream Sl_real_name, Sl_imag_name;
                Sl_real_name << "Sl_real_" << flavor1 << "_" << site1 << "_" << site2;
                Sl_imag_name << "Sl_imag_" << flavor1 << "_" << site1 << "_" << site2;
                std::vector<double> mean_real = results[Sl_real_name.str().c_str()].template mean<std::vector<double> >();
                std::vector<double> mean_imag = results[Sl_imag_name.str().c_str()].template mean<std::vector<double> >();
                for (unsigned int i_l = 0; i_l < n_legendre; ++i_l) {
                  Sl[i_l][site1][site2][flavor1] = std::complex<double>(mean_real[i_l], mean_imag[i_l]) / sign;
                }
              }
            }
          }

          //compute S(iomega_n)
          {
            Eigen::MatrixXcd Sl_vec(n_legendre, 1);
            for (int flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
              for (int site1 = 0; site1 < n_site; ++site1) {
                for (int site2 = 0; site2 < n_site; ++site2) {
                  for (int i_l = 0; i_l < n_legendre; ++i_l) {
                    Sl_vec(i_l, 0) = Sl[i_l][site1][site2][flavor1];
                  }
                  Eigen::MatrixXcd Sw_vec(Tnl * Sl_vec);
                  for (int w = 0; w < n_matsubara; ++w) {
                    Sw[w][site1][site2][flavor1] = Sw_vec(w, 0);
                  }
                }
              }
            }
          }
        }

        template<class SOLVER_TYPE>
        void compute_greens_functions(const typename alps::accumulators::result_set &results,
                                      const typename alps::params &parms, const std::string &output_file) {
          spin_t n_flavors = parms["model.flavors"];
          unsigned int n_site = parms["model.sites"];
          double beta = parms["model.beta"];
          unsigned int n_matsubara = parms["G1.n_matsubara"];
          int n_legendre = parms["G1.n_legendre"];

          std::vector<double> mean_order = results["PertOrder"].template mean<std::vector<double> >();

          std::cout << "average matrix size was: " << std::endl;
          std::ofstream matrix_size("matrix_size", std::ios::app);
          for (unsigned int i = 0; i < n_flavors; ++i) {
            std::cout << mean_order[i] << "\t";
            matrix_size << mean_order[i] << "\t";
          }
          std::cout << std::endl;
          matrix_size << std::endl;
          std::cout << "average sign was: " << results["Sign"].template mean<double>() << std::endl;

          boost::multi_array<std::complex<double>,4> Sl, Sw;

          evaluate_selfenergy_measurement_legendre<SOLVER_TYPE>(results,
                                                                beta,
                                                                n_site,
                                                                n_flavors,
                                                                n_matsubara,
                                                                n_legendre,
                                                                Sl, Sw);

          alps::hdf5::archive ar(output_file, "a");
          ar["/SigmaG_legendre"] = Sl;
          ar["/SigmaG_omega"] = Sw;
        }
    }
}
