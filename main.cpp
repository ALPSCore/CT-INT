/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Emanuel Gull <gull@phys.columbia.edu>,
 *
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include "interaction_expansion.hpp"
#include "fouriertransform.h"
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include "measurements.hpp"

typedef InteractionExpansion<real_number_solver> SOLVER_TYPE;

#ifdef ALPS_HAVE_MPI
//#include <alps/mcmpiadapter.hpp>
#include "my_mcmpiadapter.hpp"
#include "my_check_schedule.hpp"
typedef alps::mcmpiadapter<SOLVER_TYPE, alps::my_check_schedule> sim_type;
#else
typedef SOLVER_TYPE sim_type;
#endif

#undef BUILD_PYTHON_MODULE

bool stop_callback(boost::posix_time::ptime const & end_time) {
//stops the simulation if time > end_time or if signals received.
  static alps::ngs::signal signal;
  return !signal.empty() || boost::posix_time::second_clock::local_time() > end_time;
}
//void compute_greens_functions(const alps::results_type<SOLVER_TYPE>::type &results, const alps::parameters_type<SOLVER_TYPE>::type& parms, const std::string &output_file);
int global_mpi_rank;

int main(int argc, char** argv)
{
  alps::mcoptions options(argc, argv);
  if (options.valid) {
    std::string output_file = options.output_file;

#ifdef ALPS_HAVE_MPI
    //boot up MPI environment
    boost::mpi::environment env(argc, argv);
#endif

    //create ALPS parameters from hdf5 parameter file
    alps::parameters_type<SOLVER_TYPE>::type parms(alps::hdf5::archive(options.input_file, alps::hdf5::archive::READ));
    //try {
      if(options.time_limit!=0)
        throw std::invalid_argument("time limit is passed in the parameter file!");
      if(!parms.defined("MAX_TIME")) throw std::runtime_error("parameter MAX_TIME is not defined. How long do you want to run the code for? (in seconds)");

#ifndef ALPS_HAVE_MPI
      global_mpi_rank=0;
      sim_type s(parms,global_mpi_rank);
#else
      boost::mpi::communicator c;
      c.barrier();
      global_mpi_rank=c.rank();
      sim_type s(parms, c);
#endif
      //run the simulation
      s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds((int)parms["MAX_TIME"])));

      //on the master: collect MC results and store them in file, then postprocess
      if (global_mpi_rank==0){
        std::cout << " collecting ... " << std::endl;
        alps::results_type<SOLVER_TYPE>::type results = collect_results(s);
        save_results(results, parms, output_file, "/simulation/results");
        //compute the output Green's function and Fourier transform it, store in the right path
        std::cout << " collecting done rank = " << global_mpi_rank << std::endl;
        c.barrier();
        std::cout << " compute GF " << std::endl;
        compute_greens_functions<SOLVER_TYPE>(results, parms, output_file);
        std::cout << " compute GF  done" << std::endl;
#ifdef ALPS_HAVE_MPI
      } else{
        collect_results(s);
        std::cout << " collecting done rank = " << global_mpi_rank << std::endl;
        c.barrier();
      }
      c.barrier();

#else
      }
#endif
    //}
    //catch(std::exception& exc){
        //std::cerr<<exc.what()<<std::endl;
        //return -1;
    //}
    //catch(...){
        //std::cerr << "Fatal Error: Unknown Exception!\n";
        //return -2;
    //}
  }//options.valid
  return 0;
}
