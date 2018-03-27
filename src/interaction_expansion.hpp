#pragma once

#include <algorithm>
#include <fstream>
#include <cmath>

#include <alps/accumulators.hpp>
#include <alps/mc/api.hpp>
#include <alps/mc/mcbase.hpp>
#include <alps/mc/stop_callback.hpp>
#include <alps/mc/mpiadapter.hpp>

#include <boost/multi_array.hpp>
//#include <boost/timer/timer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

#include "program_options.hpp"
#include "submatrix.hpp"
#include "green_function.h"
#include "types.h"
#include "U_matrix.h"
#include "operator.hpp"
#include "legendre.h"
#include "update_statistics.h"
#include "update_manager.hpp"

namespace alps {
    namespace ctint {

/*types*/
        class c_or_cdagger;

        typedef class histogram simple_hist;

        class histogram {
        public:

            histogram(unsigned int N) : hist_(N, 0) {}

            unsigned long &operator[](unsigned int n) { return hist_[n]; }

            const unsigned long &operator[](unsigned int n) const { return hist_[n]; }

            unsigned int size() const { return hist_.size(); }

            void count(size_t k) {
              if (k < size()) {
                ++hist_[k];
              }
            }

            unsigned int max_index() const {
              unsigned int max_index = 0;
              double max = 0;
              for (unsigned int i = 0; i < hist_.size(); ++i) {
                if (max < hist_[i]) {
                  max = hist_[i];
                  max_index = i;
                }
              }
              return max_index;
            }

            unsigned int top_index() const {
              unsigned int top_index = 0;
              for (unsigned int i = 0; i < hist_.size(); ++i) {
                if (hist_[i] != 0) {
                  top_index = i;
                }
              }
              return top_index;
            }

            double max(const unsigned int index) const {
              double max = 0;
              for (unsigned int i = 0; i < index; ++i) {
                if (max < hist_[i]) {
                  max = hist_[i];
                }
              }
              return max;
            }

            double average(const unsigned int index) const {
              double average = 0;
              for (unsigned int i = 0; i < index; ++i) {
                average += hist_[i];
              }
              return average / index;
            }

            bool is_flat(const unsigned int index) const { return max(index) * 0.8 < average(index); }

            void clear() {
              for (unsigned int i = 0; i < hist_.size(); ++i) {
                hist_[i] = 0;
              }
            }

            std::valarray<double> to_valarray() {
              std::valarray<double> tmparray(hist_.size());
              for (size_t i = 0; i < hist_.size(); ++i) {
                tmparray[i] = hist_[i];
              }
              return tmparray;
            }

        private:

            std::vector<unsigned long> hist_;
        };

        typedef struct real_number_solver {
            typedef double M_TYPE;
            typedef double REAL_TYPE;
            typedef std::complex<double> COMPLEX_TYPE;
        } real_number_solver;

        typedef struct complex_number_solver {
            typedef std::complex<double> M_TYPE;
            typedef double REAL_TYPE;
            typedef std::complex<double> COMPLEX_TYPE;
        } complex_number_solver;

        /*
        template<typename T>
        class BareGreenInterpolate {
        public:
            BareGreenInterpolate(const alps::params &p);

            T operator()(const annihilator &c, const creator &cdagg) const;

            T operator()(double delta_t, int flavor, int site1, int site2) const;

            int nsite() const { return n_sites_; };

            int nflavor() const { return n_flavors_; };

            bool is_zero(int site1, int site2, int flavor, double eps) const;

        private:
            const double beta_, temp_;
            const int ntau_, n_flavors_, n_sites_;
            boost::multi_array<std::pair<T, T>, 4> AB_;
            double dbeta_;//beta/ntau
        };
        */

        template<typename T>
        std::pair<T, T> compute_weight(const general_U_matrix<T> &Uijkl, const itime_vertex_container &itime_vertices);

        class InteractionExpansionBase : public alps::mcbase {
        public:
            InteractionExpansionBase(parameters_type const &params, std::size_t seed_offset) : alps::mcbase(params,
                                                                                                                 seed_offset) {};

            virtual ~InteractionExpansionBase() {}

            virtual void update()=0;

            virtual void measure()=0;

            virtual double fraction_completed() const =0;

            // Some extension
            //virtual bool is_thermalized() const =0;

            virtual void finalize()=0;

            static parameters_type &define_parameters(parameters_type &parameters) {
              return alps::mcbase::define_parameters(parameters);
            }
        };

        template<class TYPES>
        class InteractionExpansion : public InteractionExpansionBase {
        public:

            InteractionExpansion(parameters_type const &params, std::size_t seed_offset = 42);

            ~InteractionExpansion();

            void measure();

            void update();

            double fraction_completed() const;

            void finalize();

            typedef typename TYPES::M_TYPE M_TYPE;
            typedef typename TYPES::REAL_TYPE REAL_TYPE;
            typedef typename TYPES::COMPLEX_TYPE COMPLEX_TYPE;
            typedef boost::multi_array<std::complex<double>, 4> Wk_t;
            typedef green_function<typename TYPES::COMPLEX_TYPE> matsubara_green_function_t;
            typedef green_function<typename TYPES::COMPLEX_TYPE> itime_green_function_t;

            static parameters_type &define_parameters(parameters_type &parameters) {
              InteractionExpansionBase::define_parameters(parameters);
              return define_ctint_options(parameters);
            }

        protected:

            /*functions*/
            /*io & initialization*/
            void initialize_simulation(const alps::params &parms); // called by constructor

            //void exchange_update();

            // in file io.cpp
            void print(std::ostream &os);

            // in file observables.ipp
            void measure_observables();

            void initialize_observables(void);

            void compute_Sl();

            void measure_densities();

            // in file interaction_expansion.hpp
            void sanity_check();
            //bool is_quantum_number_conserved(const itime_vertex_container& vertices);

            virtual bool is_thermalized() const {
              return step > therm_steps;
            }

            void prepare_for_measurement(); //called once after thermalization is done

            //copy of input parameters
            alps::params parms;

            /*private member variables, constant throughout the simulation*/
            //const unsigned int node;
            const unsigned int max_order;
            const spin_t n_flavors;                                //number of flavors (called 'flavors') in InteractionExpansion
            const site_t n_site;                                //number of sites
            //const itime_index_t n_tau;                        //number of imag time slices
            //const std::size_t n_legendre;
            const boost::uint64_t mc_steps;
            const boost::uint64_t therm_steps;
            const int n_ins_rem;
            const int n_shift;
            const int n_spin_flip;
            const bool force_quantum_number_conservation;
            const bool single_vertex_update_non_density_type;
            const double beta;
            //const double temperature;                        //only for performance reasons: avoid 1/beta computations where possible

            general_U_matrix<M_TYPE> Uijkl; //for any general two-body interaction

            /*heart of submatrix update*/
            typedef boost::shared_ptr<SubmatrixUpdate<M_TYPE> > WALKER_P_TYPE;
            WALKER_P_TYPE submatrix_update;

            //for measurement of Green's function
            //M is computed from A in measure_observables.
            std::vector<alps::numeric::matrix<M_TYPE> > M_flavors;

            //quantum numbers
            std::vector<std::vector<std::vector<size_t> > > groups;
            std::vector<std::vector<int> > group_map;
            std::vector<std::vector<quantum_number_t> > quantum_number_vertices;
            std::vector<int> group_dim;
            int qn_dim;

            //double vertex update
            //std::vector<std::pair<int,int> > mv_update_valid_pair;
            //boost::multi_array<bool,2> mv_update_valid_pair_flag;

            //std::vector<bool> is_density_density_type;

            //for shift update
            //std::vector<bool> shift_update_valid;

            //const unsigned int recalc_period;
            const unsigned int measurement_period;
            //const unsigned int convergence_check_period;

            /*InteractionExpansion's roundoff threshold*/
            const double almost_zero;

            bool is_thermalized_in_previous_step_;

            // Non-interacting GF
            itime_green_function_t bare_green_itime;

            //simple_hist pert_hist;
            unsigned int hist_max_index;
            simple_hist **vertex_histograms;
            unsigned int vertex_histogram_size;

            unsigned long step;
            time_t start_time;
            clock_t update_time;
            clock_t measurement_time;

            LegendreTransformer legendre_transformer;

            std::valarray<double> pert_order_hist;

            alps::mpi::communicator comm;

            //only for test
            //std::vector<typename TYPES::COMPLEX_TYPE> Wk_dynamics;
            //std::vector<typename TYPES::COMPLEX_TYPE> Sl_dynamics;
            //std::vector<double> pert_order_dynamics;

            green_function<M_TYPE> g0_intpl;

            VertexUpdateManager<M_TYPE> update_manager;

        };

/*aux functions*/
        std::ostream &operator<<(std::ostream &os, const std::vector<double> &v);

        std::ostream &operator<<(std::ostream &os, const c_or_cdagger &c);

        std::ostream &operator<<(std::ostream &os, const simple_hist &h);
    }
}


#include "interaction_expansion.ipp"
#include "measurement.ipp"

