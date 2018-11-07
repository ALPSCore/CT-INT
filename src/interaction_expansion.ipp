//#include <ctime>
//#include <boost/lexical_cast.hpp>
//#include "boost/tuple/tuple.hpp"

#include "interaction_expansion.hpp"

namespace alps {
    namespace ctint {
        template<class TYPES>
        InteractionExpansion<TYPES>::InteractionExpansion(parameters_type const &params, std::size_t seed_offset)
          : InteractionExpansionBase(params, seed_offset),
            parms(parameters),
            max_order(parms["update.max_order"]),
            n_flavors(parms["model.spins"]),
            n_site(parms["model.sites"]),
            mc_steps((boost::uint64_t) parms["total_steps"]),
            therm_steps(parms["thermalization_steps"]),
            n_ins_rem(parms["update.n_ins_rem_vertex"]),
            n_shift(parms["update.n_vertex_shift"]),
            n_spin_flip(parms["update.n_spin_flip"]),
            single_vertex_update_non_density_type(false),
            beta(parms["model.beta"]),
            Uijkl(parms),
            M_flavors(n_flavors),
            measurement_period(parms["measurement_period"].template as<int>()),
            almost_zero(1.e-16),
            is_thermalized_in_previous_step_(false),
            legendre_transformer(params["G1.n_matsubara"], params["G1.n_legendre"]),
            pert_order_hist(max_order + 1),
            comm(),
            g0_intpl(),
            update_manager(parms, Uijkl, g0_intpl, comm.rank() == 0),
            timings(4) {
          //other parameters
          step = 0;
          measurement_time = 0;
          update_time = 0;


          // Read non-interacting G(tau)
          if (params["model.G0_tau_file"] == "") {
            throw std::runtime_error("Set model.G0_tau_file!");
          }
          g0_intpl.read_itime_data(params["model.G0_tau_file"], beta, n_flavors, n_site);

          //initialize the simulation variables
          initialize_simulation(parms);

          //Some checks
          if (params["measurement_period"].template as<int>() <= 0) {
            throw std::runtime_error("measurement_period must be specified. A reasonable value is the average expansion order");
          }

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

          auto t_start = std::chrono::system_clock::now();

          std::vector<double> timing_part(2, 0.0);

          for (std::size_t i = 0; i < measurement_period; ++i) {
#ifndef NDEBUG
            std::cout << " step " << step << std::endl;
#endif
            step++;
            //std::cout << " step " << step << std::endl;

#ifndef NDEBUG
            std::cout << " step " << step << " node " << comm.rank() << " pert " << submatrix_update->pert_order() << std::endl;
#endif

            for (int i_ins_rem = 0; i_ins_rem < n_ins_rem; ++i_ins_rem) {
              update_manager.do_ins_rem_update(*submatrix_update, Uijkl, random, 1.0);
            }

            for (int i_shift = 0; i_shift < n_shift; ++i_shift) {
              update_manager.do_shift_update(*submatrix_update, Uijkl, random, false);
            }

            for (int i_spin_flip = 0; i_spin_flip < n_spin_flip; ++i_spin_flip) {
              update_manager.do_spin_flip_update(*submatrix_update, Uijkl, random);
            }

            for (spin_t flavor = 0; flavor < n_flavors; ++flavor) {
              vertex_histograms[flavor]->count(submatrix_update->invA()[flavor].creators().size());
            }

          }

          // heavy parts
          {
            auto t_start_local = std::chrono::system_clock::now();
            // I do not want to recompute the inverse matrix so frequently.
            if (step % (10 * measurement_period) == 0) {
              submatrix_update->recompute_matrix(true);
            }
            auto t_end_local = std::chrono::system_clock::now();
            timing_part[0] += std::chrono::duration_cast<std::chrono::nanoseconds>(t_end_local-t_start_local).count();
            //std::cout << "debug " << std::chrono::duration_cast<std::chrono::nanoseconds>(t_end_local-t_start_local).count() << std::endl;
          }
          {
            auto t_start_local = std::chrono::system_clock::now();
            update_manager.global_updates(submatrix_update, Uijkl, g0_intpl, random);
            auto t_end_local = std::chrono::system_clock::now();
            timing_part[1] += std::chrono::duration_cast<std::chrono::nanoseconds>(t_end_local-t_start_local).count();
          }

          if (is_thermalized() && !is_thermalized_in_previous_step_) {
            prepare_for_measurement();
          }
          is_thermalized_in_previous_step_ = is_thermalized();

          auto t_end = std::chrono::system_clock::now();

          timings[0] = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end-t_start).count()
                       - std::accumulate(timing_part.begin(), timing_part.end(), 0.0);
          timings[1] = timing_part[0];
          timings[2] = timing_part[1];
        }

        template<class TYPES>
        void InteractionExpansion<TYPES>::measure() {
          auto t_start = std::chrono::system_clock::now();

          //In the below, real physical quantities are measured.
          if (!is_thermalized()) {
            return;
          }
          measure_observables();
          update_manager.measure_observables(measurements);

          auto t_end = std::chrono::system_clock::now();
          timings[3] = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end-t_start).count();

          // from nanoseconds to milliseconds
          std::transform(timings.begin(), timings.end(), timings.begin(), [](double x){return 1E-6*x;});

          measurements["Timings"] << timings;
        }

/**
 * Finalize the Monte Carlo simulation, e.g., write some data to disk.
 */
        template<class TYPES>
        void InteractionExpansion<TYPES>::finalize() {

          /*
          std::string node_str = boost::lexical_cast<std::string>(comm.rank());

          if (parms.defined("PREFIX_DUMP_CONFIG")) {
            std::ofstream os((parms["PREFIX_DUMP_CONFIG"].template as<std::string>()
                              + std::string("-config-node") + node_str + std::string(".txt")).c_str());
            dump(os, submatrix_update->itime_vertices());
          }
           */

        }


        template<class TYPES>
        double InteractionExpansion<TYPES>::fraction_completed() const {
          if (!is_thermalized()) {
            return 0.;
          } else {
            //std::cout << "step debug " << step << " " << ((step - therm_steps) / (double) mc_steps) << std::endl;
            //std::cout << "step debug " << step << " " << step << " " <<  therm_steps << " " <<  mc_steps << std::endl;
            double fraction = (static_cast<double>(step) - static_cast<double>(therm_steps))/static_cast<double>(mc_steps);
            return fraction;
          }
        }


///do all the setup that has to be done before running the simulation.
        template<class TYPES>
        void InteractionExpansion<TYPES>::initialize_simulation(const alps::params &parms) {
          //pert_hist.clear();
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
          //std::cout << "prepare for meas" << std::endl;
          std::cout << "Rank " << comm.rank() << ": thermalization done!" << std::endl;
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
