#include <ctime>
#include <boost/lexical_cast.hpp>
#include "boost/tuple/tuple.hpp"

#include "interaction_expansion.hpp"

namespace alps {
    namespace ctint {

        /*
        template<typename T>
        BareGreenInterpolate<T>::BareGreenInterpolate(const alps::params &p) :
          beta_(p["BETA"].template as<int>()),
          temp_(1.0 / beta_),
          ntau_(p["N_TAU"].template as<int>()),
          n_flavors_(p["FLAVORS"].template as<int>()),
          n_sites_(p["SITES"].template as<int>()),
          dbeta_(beta_ / ntau_) {
          green_function<std::complex<double> > bare_green_matsubara(ntau_, n_sites_, n_flavors_),
            bare_green_itime(ntau_ + 1, n_sites_, n_flavors_);

          boost::tie(bare_green_matsubara, bare_green_itime) =
            read_bare_green_functions<std::complex<double> >(p);

          assert(ntau_ == bare_green_itime.ntime() - 1);
          AB_.resize(boost::extents[n_flavors_][n_sites_][n_sites_][ntau_ + 1]);

          for (int flavor = 0; flavor < n_flavors_; ++flavor) {
            for (int site1 = 0; site1 < n_sites_; ++site1) {
              for (int site2 = 0; site2 < n_sites_; ++site2) {
                for (int tau = 0; tau < ntau_; ++tau) {
                  const T a =
                    mycast<T>(
                      (bare_green_itime(tau + 1, site1, site2, flavor) - bare_green_itime(tau, site1, site2, flavor)) /
                      dbeta_
                    );
                  const T b = mycast<T>(bare_green_itime(tau, site1, site2, flavor));

                  AB_[flavor][site1][site2][tau] = std::make_pair(a, b);
                }
                AB_[flavor][site1][site2][ntau_] = std::make_pair(0.0, mycast<T>(
                  bare_green_itime(ntau_, site1, site2, flavor)));
              }
            }
          }
        }

        template<typename T>
        T BareGreenInterpolate<T>::operator()(const annihilator &c, const creator &cdagger) const {
          assert(c.flavor() == cdagger.flavor());

          const int flavor = c.flavor();

          const int site1 = c.s(), site2 = cdagger.s();
          double dt = c.t().time() - cdagger.t().time();
          if (dt == 0.0) {
            if (c.t().small_index() > cdagger.t().small_index()) { //G(+delta)
              return AB_[flavor][site1][site2][0].second;
            } else { //G(-delta)
              return -AB_[flavor][site1][site2][ntau_].second;
            }
          } else {
            T coeff = 1.0;
            while (dt >= beta_) {
              dt -= beta_;
              coeff *= -1.0;
            }
            while (dt < 0.0) {
              dt += beta_;
              coeff *= -1.0;
            }

            assert(dt >= 0 && dt <= beta_);
            const int time_index_1 = (int) (dt * ntau_ * temp_);
            return coeff * (AB_[flavor][site1][site2][time_index_1].first * (dt - time_index_1 * dbeta_) +
                            AB_[flavor][site1][site2][time_index_1].second);
          }
        }

//if delta_t==0, we assume delta_t = +0
        template<typename T>
        T BareGreenInterpolate<T>::operator()(double delta_t, int flavor, int site1, int site2) const {
          double dt = delta_t;
          T coeff = 1.0;
          while (dt >= beta_) {
            dt -= beta_;
            coeff *= -1.0;
          }
          while (dt < 0.0) {
            dt += beta_;
            coeff *= -1.0;
          }

          if (dt == 0.0) dt += 1E-8;

          const int time_index_1 = (int) (dt * ntau_ * temp_);
          return coeff * (AB_[flavor][site1][site2][time_index_1].first * (dt - time_index_1 * dbeta_) +
                          AB_[flavor][site1][site2][time_index_1].second);
        }

        template<typename T>
        bool BareGreenInterpolate<T>::is_zero(int site1, int site2, int flavor, double eps) const {
          return std::abs(operator()(beta_ * 1E-5, flavor, site1, site2)) < eps &&
                 std::abs(operator()(beta_ * (1 - 1E-5), flavor, site1, site2)) < eps;
        }
        */

        template<class TYPES>
        InteractionExpansion<TYPES>::InteractionExpansion(parameters_type const &params, std::size_t seed_offset)
          : InteractionExpansionBase(params, seed_offset),
            parms(parameters),
            max_order(parms["update.max_order"]),
            n_flavors(parms["model.flavors"]),
            n_site(parms["model.sites"]),
            //n_matsubara_measurements(parms["G1.NMATSUBARA_MEASUREMENTS"]),
            n_tau(parms["N_TAU"]),
            //n_tau_inv(1. / n_tau),
            //n_self(parms["NSELF"]),
            n_legendre(parms["G1.n_legendre"]),
            mc_steps((boost::uint64_t) parms["total_steps"]),
            therm_steps(parms["thermalization_steps"]),
            //max_time_in_seconds(parms["MAX_TIME"]),
            beta(parms["model.beta"]),
            temperature(1. / beta),
            Uijkl(parms),
            M_flavors(n_flavors),
            recalc_period(parms["update.recalc_period"]),
            measurement_period(parms["G1.measurement_period"].template as<int>() > 0 ? parms["G1.measurement_period"] : 500 * n_flavors * n_site),
            convergence_check_period(recalc_period),
            almost_zero(1.e-16),
            //bare_green_matsubara(n_matsubara, n_site, n_flavors),
            //bare_green_itime(n_tau + 1, n_site, n_flavors),
            pert_hist(max_order),
            legendre_transformer(params["G1.n_matsubara"], n_legendre),
            is_thermalized_in_previous_step_(false),
            n_ins_rem(parms["update.n_ins_rem_vertex"]),
            n_shift(parms["update.n_vertex_shift"]),
            n_spin_flip(parms["update.n_spin_flip"]),
            force_quantum_number_conservation(false),
            single_vertex_update_non_density_type(false),
            pert_order_hist(max_order + 1),
            comm(),
            g0_intpl(),
            update_manager(parms, Uijkl, g0_intpl, comm.rank() == 0) {

          //other parameters
          step = 0;
          start_time = time(NULL);
          measurement_time = 0;
          update_time = 0;


          // Read non-interacting G(tau)
          if (params["model.G0_tau_file"] == "") {
            throw std::runtime_error("Set model.G0_tau_file!");
          }
          g0_intpl.read_itime_data(params["model.G0_tau_file"], beta);

          //initialize the simulation variables
          initialize_simulation(parms);

          //submatrix update
          itime_vertex_container itime_vertices_init;
          if (parms.defined("PREFIX_LOAD_CONFIG")) {
            std::ifstream is((parms["PREFIX_LOAD_CONFIG"].template as<std::string>()
                              + std::string("-config-node") + boost::lexical_cast<std::string>(comm.rank()) +
                              std::string(".txt")).c_str());
            load_config<typename TYPES::M_TYPE>(is, Uijkl, itime_vertices_init);
          }

          update_manager.create_observables(measurements);


          submatrix_update = WALKER_P_TYPE(
              new SubmatrixUpdate<M_TYPE>(
                parms["update.k_ins_max"], n_flavors,
                g0_intpl, &Uijkl, beta, itime_vertices_init));

#ifndef NDEBUG
          std::cout << " step " << step << " node " << comm.rank() << " pert " << submatrix_update->pert_order() << std::endl;
#endif

          vertex_histograms = new simple_hist *[n_flavors];
          vertex_histogram_size = 100;
          for (unsigned int i = 0; i < n_flavors; ++i) {
            vertex_histograms[i] = new simple_hist(vertex_histogram_size);
          }

        }

        template<class TYPES>
        InteractionExpansion<TYPES>::~InteractionExpansion() {
          //delete submatrix_update;
        }

        template<class TYPES>
        void InteractionExpansion<TYPES>::update() {

          pert_order_hist = 0.;

          for (std::size_t i = 0; i < measurement_period; ++i) {
#ifndef NDEBUG
            std::cout << " step " << step << std::endl;
#endif
            step++;

#ifndef NDEBUG
            std::cout << " step " << step << " node " << comm.rank() << " pert " << submatrix_update->pert_order() << std::endl;
#endif

            for (int i_ins_rem = 0; i_ins_rem < n_ins_rem; ++i_ins_rem) {
              update_manager.do_ins_rem_update(*submatrix_update, Uijkl, random, 1.0);
            }

            for (int i_shift = 0; i_shift < n_shift; ++i_shift) {
              update_manager.do_shift_update(*submatrix_update, Uijkl, random, !is_thermalized());
            }

            for (int i_spin_flip = 0; i_spin_flip < n_spin_flip; ++i_spin_flip) {
              update_manager.do_spin_flip_update(*submatrix_update, Uijkl, random);
            }

            if (step % recalc_period == 0) {
              submatrix_update->recompute_matrix(true);

              update_manager.global_updates(submatrix_update, Uijkl, g0_intpl, random);
            }

            if (submatrix_update->pert_order() < max_order) {
              pert_hist[submatrix_update->pert_order()]++;
            }
            assert(submatrix_update->pert_order() < pert_order_hist.size());
            ++pert_order_hist[submatrix_update->pert_order()];

            for (spin_t flavor = 0; flavor < n_flavors; ++flavor) {
              vertex_histograms[flavor]->count(submatrix_update->invA()[flavor].creators().size());
            }

          }
        }

        template<class TYPES>
        void InteractionExpansion<TYPES>::measure() {
          //In the below, real physical quantities are measured.
          measure_observables();
          update_manager.measure_observables(measurements);
        }

/**
 * Finalize the Monte Carlo simulation, e.g., write some data to disk.
 */
        template<class TYPES>
        void InteractionExpansion<TYPES>::finalize() {

          std::string node_str = boost::lexical_cast<std::string>(comm.rank());

          if (pert_order_dynamics.size() > 0) {
            std::ofstream ofs(
              (parms["PREFIX_OUTPUT_TIME_SERIES"].template as<std::string>() + std::string("-pert_order-node") +
               node_str + std::string(".txt")).c_str());
            unsigned int n_data = pert_order_dynamics.size() / Uijkl.n_vertex_type();
            unsigned int i = 0;
            for (unsigned int i_data = 0; i_data < n_data; ++i_data) {
              ofs << i_data << " ";
              for (unsigned int i_vt = 0; i_vt < Uijkl.n_vertex_type(); ++i_vt) {
                ofs << pert_order_dynamics[i] << " ";
                ++i;
              }
              ofs << std::endl;
            }
          }

          if (Wk_dynamics.size() > 0) {
            std::ofstream ofs(
              (parms["PREFIX_OUTPUT_TIME_SERIES"].template as<std::string>() + std::string("-Wk-node") + node_str +
               std::string(".txt")).c_str());
            unsigned int n_data = Wk_dynamics.size() / (n_flavors * n_site);
            unsigned int i = 0;
            for (unsigned int i_data = 0; i_data < n_data; ++i_data) {
              ofs << i_data << " ";
              for (unsigned int flavor = 0; flavor < n_flavors; ++flavor) {
                for (unsigned int site1 = 0; site1 < n_site; ++site1) {
                  ofs << Wk_dynamics[i].real() << " " << Wk_dynamics[i].imag() << "   ";
                  ++i;
                }
              }
              ofs << std::endl;
            }
          }

          if (Sl_dynamics.size() > 0) {
            std::ofstream ofs(
              (parms["PREFIX_OUTPUT_TIME_SERIES"].template as<std::string>() + std::string("-Sl-node") + node_str +
               std::string(".txt")).c_str());
            unsigned int n_data = Sl_dynamics.size() / (n_flavors * n_site);
            unsigned int i = 0;
            for (unsigned int i_data = 0; i_data < n_data; ++i_data) {
              ofs << i_data << " ";
              for (unsigned int flavor = 0; flavor < n_flavors; ++flavor) {
                for (unsigned int site1 = 0; site1 < n_site; ++site1) {
                  ofs << Sl_dynamics[i].real() << " " << Sl_dynamics[i].imag() << "   ";
                  ++i;
                }
              }
              ofs << std::endl;
            }
          }

          if (parms.defined("PREFIX_DUMP_CONFIG")) {
            std::ofstream os((parms["PREFIX_DUMP_CONFIG"].template as<std::string>()
                              + std::string("-config-node") + node_str + std::string(".txt")).c_str());
            dump(os, submatrix_update->itime_vertices());
          }

        }


        template<class TYPES>
        double InteractionExpansion<TYPES>::fraction_completed() const {
          if (!is_thermalized()) {
            return 0.;
          } else {
            return ((step - therm_steps) / (double) mc_steps);
          }
        }


///do all the setup that has to be done before running the simulation.
        template<class TYPES>
        void InteractionExpansion<TYPES>::initialize_simulation(const alps::params &parms) {
          pert_hist.clear();
          initialize_observables();
        }

/*
template<class TYPES>
bool InteractionExpansion<TYPES>::is_quantum_number_conserved(const itime_vertex_container& vertices) {
  const int Nv = vertices.size();

  if (Nv==0)
    return true;

  std::valarray<int> qn_t(0, qn_dim), qn_max(0, qn_dim), qn_min(0, qn_dim);
  itime_vertex_container vertices_sorted(vertices);//sort vertices in decreasing order (in time)
  std::sort(vertices_sorted.begin(), vertices_sorted.end());

  for (int iv=0; iv<Nv; ++iv) {
    const vertex_definition<M_TYPE> vd = Uijkl.get_vertex(vertices_sorted[iv].type());
    vd.apply_occ_change(vertices_sorted[iv].af_state(), qn_t, qn_max, qn_min);
  }

  //check if the quantum number is conserved
  for (int i=0; i<qn_t.size(); ++i) {
    if (qn_t[i]!=0) {
      return false;
    }
  }

  return true;
}
 */


//This makes sence in the absence of a bath
/*
template<class TYPES>
bool InteractionExpansion<TYPES>::is_quantum_number_within_range(const itime_vertex_container& vertices) {
  const int Nv = vertices.size();

  if (Nv==0)
    return true;

  std::valarray<int> qn_t(0, qn_dim), qn_max(0, qn_dim), qn_min(0, qn_dim);
  itime_vertex_container vertices_sorted(vertices);//sort vertices in decreasing order (in time)
  std::sort(vertices_sorted.begin(), vertices_sorted.end());

  for (int iv=0; iv<Nv; ++iv) {
    const vertex_definition<M_TYPE> vd = Uijkl.get_vertex(vertices_sorted[iv].type());
    vd.apply_occ_change(vertices_sorted[iv].af_state(), qn_t, qn_max, qn_min);
  }

  //check if the quantum number is within range
  for (int iq=0; iq<qn_dim; ++iq) {
    //note: group_dim[iq] is zero for non-existing group
    if (qn_max[iq]-qn_min[iq]>group_dim[iq]) {
      return false;
    }
  }

  return true;
}
*/

        template<class TYPES>
        void InteractionExpansion<TYPES>::sanity_check() {
#ifndef NDEBUG
#endif
        }


        template<class TYPES>
        void InteractionExpansion<TYPES>::prepare_for_measurement() {
          update_manager.prepare_for_measurement_steps();
        }

        template<typename T, typename SPLINE_G0>
        T
        compute_weight(const general_U_matrix<T> &Uijkl, const SPLINE_G0 &spline_G0,
                       const itime_vertex_container &itime_vertices) {
          T weight_U = 1.0, weight_det = 1.0;

          const int nflavors = Uijkl.nf();
          std::vector<std::vector<annihilator> > annihilators(nflavors);
          std::vector<std::vector<creator> > creators(nflavors);
          std::vector<std::vector<T> > alpha(nflavors);
          for (int iv = 0; iv < itime_vertices.size(); ++iv) {
            const itime_vertex &v = itime_vertices[iv];
            const vertex_definition<T> &vdef = Uijkl.get_vertex(v.type());
            weight_U *= -vdef.Uval();
            for (int rank = 0; rank < v.rank(); ++rank) {
              const int flavor_rank = vdef.flavors()[rank];
              operator_time op_t(v.time(), -rank);
              creators[flavor_rank].push_back(
                creator(flavor_rank, vdef.sites()[2 * rank], op_t)
              );
              annihilators[flavor_rank].push_back(
                annihilator(flavor_rank, vdef.sites()[2 * rank + 1], op_t)
              );
              alpha[flavor_rank].push_back(vdef.get_alpha(v.af_state(), rank));
            }
          }

          alps::numeric::matrix<T> G0, work;
          for (int flavor = 0; flavor < nflavors; ++flavor) {
            const int Nv = alpha[flavor].size();
            G0.destructive_resize(Nv, Nv);
            work.destructive_resize(Nv, Nv);
            for (int j = 0; j < Nv; ++j) {
              for (int i = 0; i < Nv; ++i) {
                G0(i, j) = spline_G0(annihilators[flavor][i], creators[flavor][j]);
              }
              G0(j, j) -= alpha[flavor][j];
            }
            weight_det *= G0.safe_determinant();
          }

          return weight_U * weight_det;
        }

        template<class TYPES>
        void InteractionExpansion<TYPES>::print(std::ostream &os) {
          os
            << "***********************************************************************************************************"
            << std::endl;
          os
            << "*** InteractionExpansion solver based on ALPSCore for multi-orbital cluster impurity model              ***"
            << std::endl;
          os
            << "***     Hiroshi Shinaoka and Yusuke Nomura                                                              ***"
            << std::endl;
          os
            << "*** This code implements the interaction expansion algorithm by Rubtsov et al., JETP Letters 80, 61.    ***"
            << std::endl;
          os
            << "***                                                                                                     ***"
            << std::endl;
          os
            << "***********************************************************************************************************"
            << std::endl;
          os << parms << std::endl;
        }
    }
}
