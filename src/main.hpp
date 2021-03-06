#include "interaction_expansion.hpp"
#include "postprocess.hpp"

#ifndef ALPS_HAVE_MPI
#error ALPSCore/CT-INT requires MPI but ALPSCore was not built with MPI.
#endif

namespace alps {
    namespace ctint {

        template<class SOLVER_TYPE>
        int run_simulation(int argc, char** argv) {
          typedef alps::ctint::mcmpiadapter<SOLVER_TYPE> my_sim_type;

          //Here we construct a parameter object by parsing an ini file.
          alps::params par(argc, argv);

          char **argv_tmp = const_cast<char **>(argv);//FIXME: ugly solution
          alps::mpi::environment env(argc, argv_tmp);
          alps::mpi::communicator comm;

          my_sim_type::define_parameters(par);

          if (par.help_requested(std::cout)) {
            return 0;
          }

          if (comm.rank()==0) {
            std::cout << "Interaction-expansion QMC impurity solver based on ALPSCore" << std::endl << std::endl;
          }
          comm.barrier();

          std::cout << "Running simulation on rank " << comm.rank() << std::endl;
          my_sim_type my_sim(par, comm, 928374, par["SEED"].template as<int>());
          my_sim.run(alps::stop_callback(par["timelimit"].template as<std::size_t>()));

          // Collect the results from the simulation
          std::cout << "Rank " << comm.rank() << " has finished. Collecting results..." << std::endl;
          typename alps::results_type<my_sim_type>::type results = alps::collect_results(my_sim);

          // Print the mean and the standard deviation.
          // Only master has all the results!
          if (comm.rank()==0) {
            std::string output_file = par["outputfile"];
            alps::hdf5::archive ar(output_file, "w");
            ar["/parameters"] << par;
            ar["/simulation_raw_data/results"] << results;

            // Some post processing
            std::cout << " Postprocessing... " << std::endl;
            postprocess<SOLVER_TYPE>(results, par, output_file);
            std::cout << " done" << std::endl;
          }
          return 0;
        }
    }
}

