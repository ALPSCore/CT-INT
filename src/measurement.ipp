#include "interaction_expansion.hpp"

#include <chrono>

#include <chrono>

namespace alps {
    namespace ctint {

        inline std::vector<double> operator*(double s, const std::vector<double>& vec) {
          std::vector<double> vec2(vec);
          std::transform(vec2.begin(), vec2.end(), vec2.begin(), std::bind1st(std::multiplies<double>(), s));
          return vec2;
        }

/// Allocation of memory for the ALPS observables.
        template<class TYPES>
        void InteractionExpansion<TYPES>::initialize_observables(void) {
          measurements << SimpleRealObservable("Sign");
          measurements << SimpleRealVectorObservable("PertOrder");

          for (unsigned int flavor = 0; flavor < n_flavors; ++flavor) {
            for (unsigned int k = 0; k < n_site; k++) {
              for (unsigned int k2 = 0; k2 < n_site; k2++) {
                std::stringstream obs_name_real, obs_name_imag;
                obs_name_real << "Sl_real_" << flavor << "_" << k << "_" << k2;
                obs_name_imag << "Sl_imag_" << flavor << "_" << k << "_" << k2;
                measurements << SimpleRealVectorObservable(obs_name_real.str().c_str());
                measurements << SimpleRealVectorObservable(obs_name_imag.str().c_str());
              }
            }
          }

          measurements << SimpleRealVectorObservable("densities");
          for (int flavor = 0; flavor < n_flavors; ++flavor) {
            measurements << SimpleRealVectorObservable("densities_" + boost::lexical_cast<std::string>(flavor));
          }
          measurements << SimpleRealVectorObservable("n_i n_j");

          for (unsigned int flavor = 0; flavor < n_flavors; ++flavor) {
            for (unsigned int i = 0; i < n_site; ++i) {
              std::stringstream density_name, sz_name;
              density_name << "density_" << flavor;
              if (n_site > 1) density_name << "_" << i;
              measurements << SimpleRealObservable(density_name.str().c_str());
            }
          }

          measurements << SimpleRealVectorObservable("PerturbationOrderVertex");

          // Timings for MC update and measuremrent
          measurements << SimpleRealVectorObservable("Timings");
        }

///this function is called whenever measurements should be performed. Depending
///on the value of  measurement_method it will choose one particular
///measurement function.
        template<class TYPES>
        void InteractionExpansion<TYPES>::measure_observables() {
          //compute M from A
          submatrix_update->compute_M(M_flavors);
          const M_TYPE sign = submatrix_update->sign();

          measurements["Sign"] << alps::numeric::real(sign);
          /*
          if (parms.defined("OUTPUT_Sign") ? parms["OUTPUT_Sign"] : false) {
            std::cout << " node= " << comm.rank() << " Sign= " << sign << " pert_order= "
                      << submatrix_update->pert_order() << std::endl;
          }
          */

          measure_Sl();
          measure_densities();

          //pert_order_hist /= pert_order_hist.sum();
          //measurements["PertOrderHistogram"] << pert_order_hist;

          std::vector<double> pert_order(n_flavors);
          for (unsigned int i = 0; i < n_flavors; ++i) {
            pert_order[i] = M_flavors[i].size1();
          }
          measurements["PertOrder"] << pert_order;

          /*
          for (spin_t flavor = 0; flavor < n_flavors; ++flavor) {
            std::stringstream tmp;
            tmp << "VertexHistogram_" << flavor;
            measurements[tmp.str().c_str()] << vertex_histograms[flavor]->to_valarray();
            vertex_histograms[flavor]->clear();
          }
          */

          std::vector<double> pert_vertex(Uijkl.n_vertex_type(), 0.0);
          const itime_vertex_container &itime_vertices = submatrix_update->itime_vertices();
          for (itime_vertex_container::const_iterator it = itime_vertices.begin(); it != itime_vertices.end(); ++it) {
            assert(it->type() >= 0 && it->type() < Uijkl.n_vertex_type());
            ++pert_vertex[it->type()];
          }
          measurements["PerturbationOrderVertex"] << pert_vertex;
        }

        template<class TYPES>
        void InteractionExpansion<TYPES>::measure_Sl() {
          int n_legendre = legendre_transformer.Nl();
          int num_time_shifts = 100;

          std::vector<double> time_shifts;
          for (int i=0; i < num_time_shifts; ++i) {
            time_shifts.push_back(beta * random());
          }

          boost::multi_array<std::complex<double>, 4> Sl(boost::extents[n_flavors][n_site][n_site][n_legendre]);
          std::fill(Sl.origin(), Sl.origin() + Sl.num_elements(), 0.0);
          auto t1 = std::chrono::high_resolution_clock::now();
          compute_Sl(time_shifts, Sl);
          auto t2 = std::chrono::high_resolution_clock::now();

          boost::multi_array<std::complex<double>, 4> Sl_optimized(boost::extents[n_flavors][n_site][n_site][n_legendre]);
          std::fill(Sl_optimized.origin(), Sl_optimized.origin() + Sl_optimized.num_elements(), 0.0);
          auto t3 = std::chrono::high_resolution_clock::now();
          compute_Sl_optimized(time_shifts, Sl_optimized);
          auto t4 = std::chrono::high_resolution_clock::now();

          std::cout << "time "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t4-t3).count() << std::endl;

          auto it = Sl.origin();
          auto it2 = Sl_optimized.origin();

          for (int i=0; i<Sl.num_elements(); ++i, ++it, ++it2) {
            auto diff = std::abs(*it - *it2);
            if (diff > 1e-1) {
              throw std::runtime_error("Failed");
            }
          }

          //pass data to ALPS library
          std::vector<double> Sl_real(n_legendre, 0.0);
          std::vector<double> Sl_imag(n_legendre, 0.0);
          for (int z = 0; z < n_flavors; ++z) {
            for (int site1 = 0; site1 < n_site; ++site1) {
              for (int site2 = 0; site2 < n_site; ++site2) {
                for (int i_legendre = 0; i_legendre < n_legendre; ++i_legendre) {
                  std::complex<double> ztmp = Sl[z][site1][site2][i_legendre];
                  Sl_real[i_legendre] = ztmp.real();
                  Sl_imag[i_legendre] = ztmp.imag();
                }
                std::stringstream Sl_real_name, Sl_imag_name;
                Sl_real_name << "Sl_real_" << z << "_" << site1 << "_" << site2;
                Sl_imag_name << "Sl_imag_" << z << "_" << site1 << "_" << site2;
                measurements[Sl_real_name.str()] << Sl_real;
                measurements[Sl_imag_name.str()] << Sl_imag;
              }//site2
            }//site1
          }//z

        }

        template<class TYPES>
        void InteractionExpansion<TYPES>::compute_Sl_optimized(const std::vector<double>& time_shifts,
                                                     boost::multi_array<std::complex<double>, 4>& Sl
        ) {
          using rowmajor_mat_type = Eigen::Matrix<M_TYPE, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

          std::fill(Sl.origin(), Sl.origin() + Sl.num_elements(), 0.0);//clear the content for safety
          int n_legendre = legendre_transformer.Nl();

          M_TYPE sign = submatrix_update->sign();
          double temperature = 1.0 / beta;

          //Work arrays
          //int max_mat_size = 0;
          //for (unsigned int z = 0; z < n_flavors; ++z) {
            //max_mat_size = std::max(max_mat_size, M_flavors[z].size2());
          //}
          const std::vector<double> &sqrt_vals = legendre_transformer.get_sqrt_2l_1();

          int num_random_walk = time_shifts.size();

          std::vector<double> x_vals;
          boost::multi_array<double, 3> legendre_vals_all; //, legendre_vals_trans_all;

          for (auto z = 0; z < n_flavors; ++z) {
            int Nv = M_flavors[z].size2();

            if (Nv == 0) {
              continue;
            }

            const std::vector<annihilator> &annihilators = submatrix_update->invA()[z].annihilators();
            const std::vector<creator> &creators = submatrix_update->invA()[z].creators();

            auto t1 = std::chrono::high_resolution_clock::now();

            //interpolate G0
            boost::multi_array<M_TYPE, 3> gR(boost::extents[Nv][n_site][num_random_walk]);
            for (auto random_walk = 0; random_walk < num_random_walk; ++random_walk) {
              for (auto p = 0; p < Nv; ++p) {//annihilation operators
                double time_a = annihilators[p].t().time() + time_shifts[random_walk];
                for (auto site_B = 0; site_B < n_site; ++site_B) {
                  gR[p][site_B][random_walk] = mycast<M_TYPE>(g0_intpl(time_a, z, annihilators[p].s(), site_B));
                }
              }
            }

            auto t2 = std::chrono::high_resolution_clock::now();

            boost::multi_array<M_TYPE, 3> M_gR(boost::extents[Nv][n_site][num_random_walk]);
            Eigen::Map<rowmajor_mat_type> gR_map(gR.origin(), Nv, n_site * num_random_walk);
            Eigen::Map<rowmajor_mat_type> M_gR_map(M_gR.origin(), Nv, n_site * num_random_walk);
            M_gR_map = M_flavors[z].block() * gR_map;

            auto t3 = std::chrono::high_resolution_clock::now();

            // Compute values of Legendre polynomials
            x_vals.resize(0);//num_random_walk * Nv);
            legendre_vals_all.resize(boost::extents[n_legendre][num_random_walk][Nv]);
            for (auto random_walk = 0; random_walk < num_random_walk; ++random_walk) {
              for (auto q = 0; q < Nv; ++q) {//creation operators
                double tmp = creators[q].t().time() + time_shifts[random_walk];
                double time_c_shifted = tmp < beta ? tmp : tmp - beta;
                x_vals.push_back(2 * time_c_shifted * temperature - 1.0);
              }
            }
            auto ref = boost::multi_array_ref<double,2>(legendre_vals_all.origin(), boost::extents[n_legendre][num_random_walk * Nv]);
            legendre_transformer.compute_legendre(x_vals, ref);//P_l[x(tau_q)]
            auto t4 = std::chrono::high_resolution_clock::now();

            // Compute Sl
            for (auto random_walk = 0; random_walk < num_random_walk; ++random_walk) {
              for (auto q = 0; q < Nv; ++q) {//creation operators
                auto site_c = creators[q].s();
                double tmp = creators[q].t().time() + time_shifts[random_walk];
                std::complex<double> coeff = tmp < beta ? 1 : -1;

                coeff *= sign / (1. * num_random_walk);

                for (unsigned int site_B = 0; site_B < n_site; ++site_B) {
                  for (unsigned int i_legendre = 0; i_legendre < n_legendre; ++i_legendre) {
                    Sl[z][site_c][site_B][i_legendre] +=
                      coeff * sqrt_vals[i_legendre] * legendre_vals_all[i_legendre][random_walk][q] * M_gR[q][site_B][random_walk];
                  }
                }
              }
            }
            auto t5 = std::chrono::high_resolution_clock::now();

            std::cout << "t2 - t1 " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count() << std::endl;
            std::cout << "t3 - t2 " << std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count() << std::endl;
            std::cout << "t4 - t3 " << std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count() << std::endl;
            std::cout << "t5 - t4 " << std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count() << std::endl;

          }
        }

        template<class TYPES>
        void InteractionExpansion<TYPES>::compute_Sl(const std::vector<double>& time_shifts,
                                                     boost::multi_array<std::complex<double>, 4>& Sl
          ) {
          std::fill(Sl.origin(), Sl.origin() + Sl.num_elements(), 0.0);//clear the content for safety
          int n_legendre = legendre_transformer.Nl();

          const M_TYPE sign = submatrix_update->sign();
          const double temperature = 1.0 / beta;

          //Work arrays
          int max_mat_size = 0;
          for (unsigned int z = 0; z < n_flavors; ++z) {
            max_mat_size = std::max(max_mat_size, M_flavors[z].size2());
          }
          const std::vector<double> &sqrt_vals = legendre_transformer.get_sqrt_2l_1();

          int num_random_walk = time_shifts.size();

          std::vector<double> x_vals;
          boost::multi_array<double, 2> legendre_vals_all; //, legendre_vals_trans_all;

          alps::numeric::matrix<M_TYPE> gR(max_mat_size, n_site), M_gR(max_mat_size, n_site);

          for (unsigned int z = 0; z < n_flavors; ++z) {
            int Nv = M_flavors[z].size2();

            if (Nv == 0) {
              continue;
            }
            gR.destructive_resize(Nv, n_site);
            M_gR.destructive_resize(Nv, n_site);

            x_vals.resize(Nv);
            legendre_vals_all.resize(boost::extents[n_legendre][Nv]);

            const std::vector<annihilator> &annihilators = submatrix_update->invA()[z].annihilators();
            const std::vector<creator> &creators = submatrix_update->invA()[z].creators();

            //shift times of operators by time_shift
            for (std::size_t random_walk = 0; random_walk < num_random_walk; ++random_walk) {

              auto t1 = std::chrono::high_resolution_clock::now();

              double time_shift = time_shifts[random_walk];

              for (unsigned int p = 0; p < Nv; ++p) {//annihilation operators
                const double time_a = annihilators[p].t().time() + time_shift;

                //interpolate G0
                for (unsigned int site_B = 0; site_B < n_site; ++site_B) {
                  gR(p, site_B) = mycast<M_TYPE>(g0_intpl(time_a, z, annihilators[p].s(), site_B));
                }
              }

              auto t2 = std::chrono::high_resolution_clock::now();
              gemm(M_flavors[z], gR, M_gR);

              auto t3 = std::chrono::high_resolution_clock::now();
              //compute legendre coefficients
              for (unsigned int q = 0; q < Nv; ++q) {//creation operators
                const double tmp = creators[q].t().time() + time_shift;
                const double time_c_shifted = tmp < beta ? tmp : tmp - beta;
                x_vals[q] = 2 * time_c_shifted * temperature - 1.0;
              }
              legendre_transformer.compute_legendre(x_vals, legendre_vals_all);//P_l[x(tau_q)]

              auto t4 = std::chrono::high_resolution_clock::now();
              for (unsigned int q = 0; q < Nv; ++q) {//creation operators
                const unsigned int site_c = creators[q].s();
                const double tmp = creators[q].t().time() + time_shift;
                std::complex<double> coeff = tmp < beta ? 1 : -1;

                coeff *= sign / (1. *num_random_walk);

                for (unsigned int site_B = 0; site_B < n_site; ++site_B) {
                  for (unsigned int i_legendre = 0; i_legendre < n_legendre; ++i_legendre) {
                    Sl[z][site_c][site_B][i_legendre] +=
                      coeff * sqrt_vals[i_legendre] * legendre_vals_all[i_legendre][q] * M_gR(q, site_B);
                  }
                }
              }
              auto t5 = std::chrono::high_resolution_clock::now();

              //std::cout << "t2 - t1 " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count() << std::endl;
              //std::cout << "t3 - t2 " << std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count() << std::endl;
              //std::cout << "t4 - t3 " << std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count() << std::endl;
              //std::cout << "t5 - t4 " << std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count() << std::endl;
            }//random_walk
          }
        }


        template<class TYPES>
        void InteractionExpansion<TYPES>::measure_densities() {
          const M_TYPE sign = submatrix_update->sign();

          std::vector<std::vector<double> > dens(n_flavors);
          for (unsigned int z = 0; z < n_flavors; ++z) {
            dens[z].resize(n_site);
            memset(&(dens[z][0]), 0., sizeof(double) * (n_site));
          }
          double tau = beta * random();
          double sign_real = mycast<double>(sign);
          for (unsigned int z = 0; z < n_flavors; ++z) {
            const size_t Nv = M_flavors[z].size2();
            const std::vector<annihilator> &annihilators = submatrix_update->invA()[z].annihilators();
            const std::vector<creator> &creators = submatrix_update->invA()[z].creators();

            using eigen_vector_t = Eigen::Matrix<M_TYPE, Eigen::Dynamic, 1>;
            eigen_vector_t g0_tauj(Nv), M_g0_tauj(Nv), g0_taui(Nv);

            for (unsigned int s = 0; s < n_site; ++s) {
              for (unsigned int j = 0; j < Nv; ++j)
                g0_tauj[j] = mycast<M_TYPE>(g0_intpl(annihilators[j].t().time() - tau, z, annihilators[j].s(),
                                                     s));//CHECK THE TREATMENT OF EQUAL-TIME Green's function
              for (unsigned int i = 0; i < Nv; ++i)
                g0_taui[i] = mycast<M_TYPE>(g0_intpl(tau - creators[i].t().time(), z, s, creators[i].s()));
              if (M_flavors[z].size2() > 0) {
                M_g0_tauj = M_flavors[z].block() * g0_tauj;
              }
              dens[z][s] += mycast<double>(g0_intpl(-beta * 1E-10, z, s, s));//tau=-0
              for (unsigned int j = 0; j < Nv; ++j) {
                dens[z][s] -= mycast<double>(g0_taui[j] * M_g0_tauj[j]);
              }
            }
          }
          std::vector<double> densities(n_flavors, 0.0);
          for (unsigned int z = 0; z < n_flavors; ++z) {
            std::vector<double> densmeas(n_site);
            for (unsigned int i = 0; i < n_site; ++i) {
              densities[z] += dens[z][i];
              densmeas[i] = dens[z][i];
            }
            {
              std::vector<double> signed_densmeas(densmeas);
              for (std::size_t i = 0; i < n_site; ++i) {
                signed_densmeas[i] *= sign_real;
              }
            }
            measurements["densities_" + boost::lexical_cast<std::string>(z)] << sign_real * densmeas;
            densities[z] /= n_site;
            densities[z] = densities[z];
          }
          measurements["densities"] << sign_real * densities;

          {
            std::vector<double> ninj(n_site * n_site * n_flavors * n_flavors);
            int pos = 0;
            for (unsigned int flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
              for (unsigned int i = 0; i < n_site; ++i) {
                for (unsigned int flavor2 = 0; flavor2 < n_flavors; ++flavor2) {
                  for (unsigned int j = 0; j < n_site; ++j) {
                    ninj[pos] = (dens[flavor1][i]) * (dens[flavor2][j]);
                    ++pos;
                  }
                }
              }
            }
            measurements["n_i n_j"] << sign_real * ninj;
          }
        }
    }
}
