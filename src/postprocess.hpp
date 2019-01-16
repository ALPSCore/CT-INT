//
// Created by H. Shinaoka on 2015/12/23.
//

#pragma once

#include <complex>

#include <alps/accumulators.hpp>

#include "legendre.h"
#include "hdf5/boost_any.hpp"

namespace alps {
    namespace ctint {

        template<class SOLVER_TYPE>
        void evaluate_selfenergy_measurement_legendre(const typename alps::accumulators::result_set &results,
                                                      const typename alps::params &parms,
                                                      boost::multi_array<std::complex<double>,4>& Sl,
                                                      boost::multi_array<std::complex<double>,4>& Sw
        ) {
          std::cout << "evaluating self energy measurement: lengendre, real space" << std::endl;
          double sign = results["Sign"].template mean<double>();

          int n_site = parms["model.sites"];
          int n_flavors = parms["model.spins"];
          int n_matsubara = parms["G1.n_matsubara"];
          int n_legendre = parms["G1.n_legendre"];

          //Legendre expansion utils
          LegendreTransformer legendre_transformer(n_matsubara, n_legendre);
          auto &Tnl = legendre_transformer.Tnl();

          Sl.resize(boost::extents[n_legendre][n_site][n_site][n_flavors]);
          Sw.resize(boost::extents[n_matsubara][n_site][n_site][n_flavors]);

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
        void evaluate_selfenergy_measurement_binning(const typename alps::accumulators::result_set &results,
                                                     const typename alps::params &parms,
                                                     boost::multi_array<std::complex<double>,4>& Stau,
                                                     boost::multi_array<std::complex<double>,4>& Sl,
                                                     boost::multi_array<std::complex<double>,4>& Sw
        ) {
          std::cout << "evaluating self energy measurement: binning, real space" << std::endl;
          double sign = results["Sign"].template mean<double>();

          double beta = parms["model.beta"];
          int n_site = parms["model.sites"];
          int n_flavors = parms["model.spins"];
          int n_matsubara = parms["G1.n_matsubara"];
          int n_legendre = parms["G1.n_legendre"];
          int n_bin = parms["G1.n_bin"];

          Stau.resize(boost::extents[n_bin][n_site][n_site][n_flavors]);
          Sl.resize(boost::extents[n_legendre][n_site][n_site][n_flavors]);
          std::fill(Sl.origin(), Sl.origin() + Sl.num_elements(), 0.0);
          Sw.resize(boost::extents[n_matsubara][n_site][n_site][n_flavors]);

          //load Stau
          for (auto flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
            for (auto site1 = 0; site1 < n_site; ++site1) {
              for (auto site2 = 0; site2 < n_site; ++site2) {
                std::stringstream Stau_real_name, Stau_imag_name;
                Stau_real_name << "Stau_real_" << flavor1 << "_" << site1 << "_" << site2;
                Stau_imag_name << "Stau_imag_" << flavor1 << "_" << site1 << "_" << site2;
                std::vector<double> mean_real = results[Stau_real_name.str().c_str()].template mean<std::vector<double> >();
                std::vector<double> mean_imag = results[Stau_imag_name.str().c_str()].template mean<std::vector<double> >();
                for (auto ib = 0; ib < n_bin; ++ib) {
                  Stau[ib][site1][site2][flavor1] = std::complex<double>(mean_real[ib], mean_imag[ib]) / sign;
                }
              }
            }
          }

          //Transform to Legendre basis
          LegendreTransformer legendre_transformer(n_matsubara, n_legendre);
          const std::vector<double> &sqrt_vals = legendre_transformer.get_sqrt_2l_1();
          std::vector<double> x_vals(n_bin);
          auto dtau = beta/n_bin;
          auto dx = 2.0/n_bin;
          for (auto ib=0; ib < n_bin; ++ib) {
            x_vals[ib] = (ib + 0.5) * dx - 1;
          }
          boost::multi_array<double,2> Plx(boost::extents[n_legendre][n_bin]);
          legendre_transformer.compute_legendre(x_vals, Plx);
          //for (int i_l = 0; i_l < n_legendre; ++i_l) {
            //for (auto ib = 0; ib < n_bin; ++ib) {
              //std::cout << "debug " << i_l << " " << ib << " " << Plx[i_l][ib] << std::endl;
            //}
          //}

          for (auto flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
            for (auto site1 = 0; site1 < n_site; ++site1) {
              for (auto site2 = 0; site2 < n_site; ++site2) {
                for (int i_l = 0; i_l < n_legendre; ++i_l) {
                  for (auto ib = 0; ib < n_bin; ++ib) {
                    Sl[i_l][site1][site2][flavor1] += sqrt_vals[i_l] * Plx[i_l][ib] * Stau[ib][site1][site2][flavor1] * dtau;
                  }
                }
              }
            }
          }

          //compute S(iomega_n)
          {
            auto &Tnl = legendre_transformer.Tnl();
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
        void postprocess_densities(const typename alps::accumulators::result_set &results,
                                                      const typename alps::params &parms, alps::hdf5::archive& ar) {
          std::cout << "evaluating densities" << std::endl;
          double sign = results["Sign"].template mean<double>();

          int n_site = parms["model.sites"];
          int n_flavors = parms["model.spins"];

          auto divide_by_sign =  [&](double x){return x/sign;};
          std::vector<double> densities = results["densities"].template mean<std::vector<double> >();
          std::transform(densities.begin(), densities.end(), densities.begin(), divide_by_sign);
          ar["/densities"] = densities;

          std::vector<double> ni_nj_flatten = results["n_i n_j"].template mean<std::vector<double> >();
          boost::multi_array<double,4> ni_nj(boost::extents[n_flavors][n_site][n_flavors][n_site]);
          std::transform(ni_nj_flatten.begin(), ni_nj_flatten.end(), ni_nj.origin(), divide_by_sign);
          ar["/ni_nj"] = ni_nj;
        }

        template<class SOLVER_TYPE>
        void postprocess(const typename alps::accumulators::result_set &results,
                                      const typename alps::params &parms, const std::string &output_file) {
          spin_t n_flavors = parms["model.spins"];

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

          alps::hdf5::archive ar(output_file, "a");

          /*  Single-particle Green's function */
          //boost::multi_array<std::complex<double>,4> Stau, Sl, Sw;
          boost::multi_array<std::complex<double>,4> Stau, Sl_binning, Sw_binning;
          //evaluate_selfenergy_measurement_legendre<SOLVER_TYPE>(results, parms, Sl, Sw);
          evaluate_selfenergy_measurement_binning<SOLVER_TYPE>(results, parms, Stau, Sl_binning, Sw_binning);
          ar["Sign"] = results["Sign"].template mean<double>();
          //ar["/SigmaG_legendre"] = Sl;
          //ar["/SigmaG_omega"] = Sw;
          ar["/SigmaG_legendre"] = Sl_binning;
          ar["/SigmaG_omega"] = Sw_binning;
          ar["/SigmaG_tau"] = Stau;


          /* Density and density correlations */
          postprocess_densities<SOLVER_TYPE>(results, parms, ar);

          /* Timings */
          std::vector<double> timings = results["Timings"].template mean<std::vector<double> >();
          std::cout << std::endl << "#### Timing analysis ####" << std::endl;
          std::cout << "For measurement_period = " << parms["measurement_period"].template as<int>()
                    << " steps, each part took " << std::endl
                    << " Monte Carlo update: "  << timings[0] << " ms" << std::endl
                    << " Recompute inverse matrix: "  << timings[1] << " ms" << std::endl
                    << " Global update: "  << timings[2] << " ms" << std::endl
                    << " Measurement: "  << timings[3] << " ms" << std::endl
                    << std::endl;
          //std::cout << "If the latter dominates, please increase the value of measurement_period." << std::endl << std::endl;
        }
    }
}
