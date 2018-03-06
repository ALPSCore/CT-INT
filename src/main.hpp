#include "interaction_expansion.hpp"
#include "fouriertransform.h"
#include "measurements.hpp"

#ifndef ALPS_HAVE_MPI
#error ALPSCore/CT-INT requires MPI but ALPSCore was not built with MPI.
#endif

namespace alps {
    namespace ctint {

        template<class SOLVER_TYPE>
        int run_simulation(int argc, char** argv) {
          typedef alps::mcmpiadapter<SOLVER_TYPE> my_sim_type;

          //Here we construct a parameter object by parsing an ini file.
          alps::params par(argc, argv);

          char **argv_tmp = const_cast<char **>(argv);//FIXME: ugly solution
          alps::mpi::environment env(argc, argv_tmp);
          alps::mpi::communicator comm;

          my_sim_type::define_parameters(par);
          if (par.help_requested(std::cout)) {
            return 0;
          }

          std::cout << "Running simulation on rank " << comm.rank() << std::endl;
          my_sim_type my_sim(par,comm);
          my_sim.run(alps::stop_callback(5));

          // Collect the results from the simulation
          std::cout << "Rank " << comm.rank() << " has finished. Collecting results..." << std::endl;
          typename alps::results_type<my_sim_type>::type results = alps::collect_results(my_sim);

          // Print the mean and the standard deviation.
          // Only master has all the results!
          /*
          if (comm.rank()==0) {
            std::string output_file = parameters["outputfile"];
            alps::hdf5::archive ar(output_file, "w");
            ar["/simulation/results"] << results;

            // Some post processing
            std::cout << " computing GF " << std::endl;
            compute_greens_functions<SOLVER_TYPE>(results, parms, output_file);
            std::cout << " compute GF done" << std::endl;
          }
          */
          return 0;
        }
    }
}

