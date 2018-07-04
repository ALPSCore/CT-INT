#pragma once

#include <alps/params.hpp>
#include <alps/utilities/fs/remove_extensions.hpp>

namespace alps {
    namespace ctint {
        inline alps::params& define_ctint_options(alps::params &parms) {
          parms.description("Continous-time interaction expansion impurity solver");
          parms.define<long>("total_steps", 0, "Number of Monte Carlo sweeps");
          parms.define<std::size_t>("timelimit", 0, "Total simulation time (in units of second). 0 means \"indefinitely\"");
          parms.define<long>("thermalization_steps", 0, "Number of thermalization steps");
          parms.define<int>("measurement_period", -1, "Interval between measurements");
          parms.define<std::string>("outputfile", alps::fs::remove_extensions(origin_name(parms)) + ".out.h5", "name of the output file");

          //model
          parms.define<int>("model.sites", "Number of sites");
          parms.define<int>("model.spins", "Number of spins per site");//internally, model.spins is denoted by "flavors".
          parms.define<std::string>("model.U_matrix_file", "Text file containing a list of interaction terms");
          parms.define<double>("model.U", "onsite U");
          parms.define<std::string>("model.G0_tau_file", "", "Text file containing non-interacting Green's function");
          parms.define<double>("model.beta", "Inverse temperature");

          //update
          parms.define<int>("update.max_order", 10240, "Max perturbation order");
          parms.define<int>("update.n_ins_rem_vertex", 1, "????? ");
          parms.define<int>("update.n_vertex_shift", 1, "How many vertex shift updates are performed at each MC step");
          parms.define<double>("update.vertex_shift_step_size", 0.1, "Step size for shift updates in units of beta");
          parms.define<int>("update.n_spin_flip", 1, "How many spin flip updates are performed at each MC step");
          parms.define<int>("update.k_ins_max", 100, "Batch size for submatrix update: k^ins_max in PRB 89, 195146 (2014)");
          parms.define<int>("update.n_multi_vertex_update", 1, "????? ");
          parms.define<int>("update.n_tau_statistics", 100, "Number of tau points for statistics");

          //Measurement
          parms.define<int>("G1.n_legendre", 200, "Number of Legendre polynomials");
          parms.define<int>("G1.n_matsubara", 1024, "Number of Matsubara frequencies");

          //parms.define<int>("MAX_TIME", 86400, "Max simulation time in units of second");

          //parms.define<bool>("FORCE_QUANTUM_NUMBER_CONSERVATION", false, "Will be removed.");
          //parms.define<bool>(single_vertex_update_non_density_type(parms.defined("SINGLE_VERTEX_UPDATE_FOR_NON_DENSITY_TYPE") ? parms["SINGLE_VERTEX_UPDATE_FOR_NON_DENSITY_TYPE"] : true),

          return parms;
        }
    }
}
