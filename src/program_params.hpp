#pragma once
#include <alps/params.hpp>

namespace alps {
    namespace ctint {
        inline alps::params& define_ctint_parameters(alps::params &parms) {
          parms.define<int>("MAX_ORDER", 10240, "Max perturbation order");
          parms.define<int>("FLAVORS", "Number of flavors per site");
          parms.define<int>("SITES", "Number of sites");
          parms.define<int>("NMATSUBARA", "Number of Matsubara frequencies");
          parms.define<int>("NMATSUBARA_MEASUREMENTS", 0, "??");
          parms.define<int>("N_TAU", "Number of imaginary times");
          //parms.define<int>("NSELF", "UNKNOWN?");
          parms.define<int>("N_LEGENDRE", 0, "Number of Legendre polynomials");
          parms.define<long>("SWEEPS", 0, "Number of Monte Carlo sweeps");
          parms.define<long>("THERMALIZATION", 0, "Number of thermalization steps");
          parms.define<int>("MAX_TIME", 86400, "Max simulation time in units of second");
          parms.define<double>("BETA", "Inverse temperature");
          parms.define<int>("RECALC_PERIOD", 5000, "Interval for recomputing determinat matrix from scratch");
          parms.define<int>("N_INS_REM_VERTEX", 1, "????? ");
          parms.define<int>("N_VERTEX_SHIFT", 1, "How often vertex shift is performed");
          parms.define<int>("N_SPIN_FLIP", 1, "How often spin flip is performed");
          parms.define<bool>("FORCE_QUANTUM_NUMBER_CONSERVATION", false, "AA");
          parms.define<int>("K_INS_MAX", 32, "Batch size for submatrix update");
          parms.define<int>("N_MULTI_VERTEX_UPDATE", 1, "????? ");
          parms.define<int>("N_TAU_UPDATE_STATISTICS", 100, "????? ");
          //parms.define<bool>(single_vertex_update_non_density_type(parms.defined("SINGLE_VERTEX_UPDATE_FOR_NON_DENSITY_TYPE") ? parms["SINGLE_VERTEX_UPDATE_FOR_NON_DENSITY_TYPE"] : true),

          return parms;
        }
    }
}
