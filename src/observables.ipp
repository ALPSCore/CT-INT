#include "interaction_expansion.hpp"

namespace alps {
    namespace ctint {

#include "accumulators.hpp"
#include <complex>

///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables. It is also called at the start of every DMFT
///iteration.
        template<class TYPES>
        void InteractionExpansion<TYPES>::initialize_observables(void) {
          if (measurements.has("Sign")) {
            measurements.clear();
          }
          measurements << SimpleRealObservable("Sign");
          measurements << SimpleRealVectorObservable("PertOrder");
          if (n_matsubara_measurements > 0) {
            for (unsigned int flavor = 0; flavor < n_flavors; ++flavor) {
              for (unsigned int k = 0; k < n_site; k++) {
                for (unsigned int k2 = 0; k2 < n_site; k2++) {
                  std::stringstream obs_name_real, obs_name_imag;
                  obs_name_real << "Wk_real_" << flavor << "_" << k << "_" << k2;
                  obs_name_imag << "Wk_imag_" << flavor << "_" << k << "_" << k2;
                  measurements << SimpleRealVectorObservable(obs_name_real.str().c_str());
                  measurements << SimpleRealVectorObservable(obs_name_imag.str().c_str());
                }
              }
            }
          }

          if (n_legendre > 0) {
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
          }

          measurements << SimpleRealVectorObservable("densities");
          for (int flavor = 0; flavor < n_flavors; ++flavor) {
            measurements << SimpleRealVectorObservable("densities_" + boost::lexical_cast<std::string>(flavor));
          }
          measurements << SimpleRealObservable("density_correlation");
          measurements << SimpleRealVectorObservable("n_i n_j");

          for (unsigned int flavor = 0; flavor < n_flavors; ++flavor) {
            for (unsigned int i = 0; i < n_site; ++i) {
              std::stringstream density_name, sz_name;
              density_name << "density_" << flavor;
              if (n_site > 1) density_name << "_" << i;
              measurements << SimpleRealObservable(density_name.str().c_str());
            }
          }
          for (unsigned int i = 0; i < n_site; ++i) {
            std::stringstream sz_name, sz2_name, sz0_szj_name;
            sz_name << "Sz_" << i;
            sz2_name << "Sz2_" << i;
            sz0_szj_name << "Sz0_Sz" << i;
          }
          measurements << SimpleRealVectorObservable("MeasurementTimeMsec");
          measurements << SimpleRealVectorObservable("UpdateTimeMsec");
          measurements << SimpleRealVectorObservable("UpdateTimeMsecAllWalkers");
          measurements << SimpleRealObservable("RecomputeTime");
          for (spin_t flavor = 0; flavor < n_flavors; ++flavor) {
            std::stringstream tmp;
            tmp << "VertexHistogram_" << flavor;
            measurements << SimpleRealVectorObservable(tmp.str().c_str());
          }

          measurements << SimpleRealVectorObservable("PerturbationOrderVertex");
          measurements << SimpleRealVectorObservable("PertOrderHistogram");
          measurements << SimpleRealVectorObservable("ACCEPTANCE_RATE_EXCHANGE");
        }

///this function is called whenever measurements should be performed. Depending
///on the value of  measurement_method it will choose one particular
///measurement function.
        template<class TYPES>
        void InteractionExpansion<TYPES>::measure_observables() {
          //compute M from A
          submatrix_update->compute_M(M_flavors);
          const M_TYPE sign = submatrix_update->sign();

          measurements["Sign"] << mycast<REAL_TYPE>(sign);
          if (parms.defined("OUTPUT_Sign") ? parms["OUTPUT_Sign"] : false) {
            std::cout << " node= " << comm.rank() << " Sign= " << sign << " pert_order= "
                      << submatrix_update->pert_order() << std::endl;
          }

          pert_order_hist /= pert_order_hist.sum();
          measurements["PertOrderHistogram"] << pert_order_hist;

          //const double t1 = timer.elapsed().wall*1E-6;
          if (n_matsubara_measurements > 0) {
            compute_W_matsubara();
          }
          //const double t2 = timer.elapsed().wall*1E-6;
          if (n_legendre > 0) {
            compute_Sl();
          }

          std::valarray<double> pert_order(n_flavors);
          for (unsigned int i = 0; i < n_flavors; ++i) {
            pert_order[i] = M_flavors[i].size1();
          }
          measurements["PertOrder"] << pert_order;

          for (spin_t flavor = 0; flavor < n_flavors; ++flavor) {
            std::stringstream tmp;
            tmp << "VertexHistogram_" << flavor;
            measurements[tmp.str().c_str()] << vertex_histograms[flavor]->to_valarray();
            vertex_histograms[flavor]->clear();
          }

          std::valarray<double> pert_vertex(Uijkl.n_vertex_type());
          const itime_vertex_container &itime_vertices = submatrix_update->itime_vertices();
          for (itime_vertex_container::const_iterator it = itime_vertices.begin(); it != itime_vertices.end(); ++it) {
            assert(it->type() >= 0 && it->type() < Uijkl.n_vertex_type());
            ++pert_vertex[it->type()];
          }
          measurements["PerturbationOrderVertex"] << pert_vertex;
        }
    }
}
