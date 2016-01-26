/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
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
#include <complex>
#include <alps/alea.h>
#include <alps/alea/simpleobseval.h>
#include <alps/scheduler/montecarlo.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>

///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables. It is also called at the start of every DMFT
///iteration.
template<class TYPES>
void InteractionExpansion<TYPES>::initialize_observables(void)
{
  if(measurements.has("Sign")){
    measurements.clear();
  }
#ifdef ALPS_NGS_USE_NEW_ALEA
  measurements << alps::accumulator::RealObservable("Sign");
  measurements << alps::accumulator::RealVectorObservable("PertOrder");
#else
  measurements << alps::ngs::RealObservable("Sign");
  measurements << alps::ngs::RealVectorObservable("PertOrder");
#endif
  if(n_matsubara_measurements>0) {
    for(unsigned int flavor=0;flavor<n_flavors;++flavor){
      for (unsigned int k = 0; k < n_site; k++) {
        for (unsigned int k2 = 0; k2 < n_site; k2++) {
          std::stringstream obs_name_real, obs_name_imag;
          obs_name_real << "Wk_real_" << flavor << "_" << k << "_" << k2;
          obs_name_imag << "Wk_imag_" << flavor << "_" << k << "_" << k2;
#ifndef ALPS_NGS_USE_NEW_ALEA
          measurements << alps::ngs::RealVectorObservable(obs_name_real.str().c_str());
          measurements << alps::ngs::RealVectorObservable(obs_name_imag.str().c_str());
#else
          throw std::runtime_error("alps::ngs::SignedRealVectorObservable is not implemented");
#endif //ALPS_NGS_USE_NEW_ALEA
        }
      }
    }
  }

  if (n_legendre>0) {
    for(unsigned int flavor=0;flavor<n_flavors;++flavor){
      for (unsigned int k = 0; k < n_site; k++) {
        for (unsigned int k2 = 0; k2 < n_site; k2++) {
          std::stringstream obs_name_real, obs_name_imag;
          obs_name_real << "Sl_real_" << flavor << "_" << k << "_" << k2;
          obs_name_imag << "Sl_imag_" << flavor << "_" << k << "_" << k2;
#ifndef ALPS_NGS_USE_NEW_ALEA
          measurements << alps::ngs::RealVectorObservable(obs_name_real.str().c_str());
          measurements << alps::ngs::RealVectorObservable(obs_name_imag.str().c_str());
#else
          throw std::runtime_error("alps::ngs::SignedRealVectorObservable is not implemented");
#endif //ALPS_NGS_USE_NEW_ALEA
        }
      }
    }
  }

#ifndef ALPS_NGS_USE_NEW_ALEA
  measurements << alps::ngs::RealVectorObservable("densities");
  for(unsigned int flavor=0;flavor<n_flavors;++flavor)
    measurements << alps::ngs::RealVectorObservable("densities_"+boost::lexical_cast<std::string>(flavor));
  measurements << alps::ngs::RealObservable("density_correlation");
  measurements << alps::ngs::RealVectorObservable("n_i n_j");
#else
  throw std::runtime_error("alps::ngs::SignedRealVectorObservable is not implemented");
#endif //ALPS_NGS_USE_NEW_ALEA
  for(unsigned int flavor=0;flavor<n_flavors;++flavor){
    for(unsigned int i=0;i<n_site;++i){
      std::stringstream density_name, sz_name;
      density_name<<"density_"<<flavor;
      if (n_site>1) density_name<<"_"<<i;
#ifndef ALPS_NGS_USE_NEW_ALEA
      measurements << alps::ngs::RealObservable(density_name.str().c_str());
#else
  throw std::runtime_error("alps::ngs::SignedRealVectorObservable is not implemented");
#endif //ALPS_NGS_USE_NEW_ALEA
    }
  }
  for(unsigned int i=0;i<n_site;++i){
    std::stringstream sz_name, sz2_name, sz0_szj_name;
    sz_name<<"Sz_"<<i;
    sz2_name<<"Sz2_"<<i;
    sz0_szj_name<<"Sz0_Sz"<<i;
// #ifndef ALPS_NGS_USE_NEW_ALEA
//     measurements << alps::ngs::SignedRealObservable(sz_name.str().c_str());
//     measurements << alps::ngs::SignedRealObservable(sz2_name.str().c_str());
//     measurements << alps::ngs::SignedRealObservable(sz0_szj_name.str().c_str());
// #else
   //throw std::runtime_error("alps::ngs::SignedRealVectorObservable is not implemented");
// #endif //ALPS_NGS_USE_NEW_ALEA
  }
  //acceptance probabilities
#ifdef ALPS_NGS_USE_NEW_ALEA
  measurements << alps::accumulator::RealObservable("VertexInsertion");
  measurements << alps::accumulator::RealObservable("VertexRemoval");
  measurements << alps::accumulator::RealObservable("MeasurementTime");
  measurements << alps::accumulator::RealObservable("UpdateTime");
  measurements << alps::accumulator::RealObservable("RecomputeTime");
  measurements << alps::accumulator::SimpleRealObservable("VertexHistogram");
#else
  for (int iv=0; iv<n_multi_vertex_update; ++iv)  {
      measurements << alps::ngs::RealVectorObservable("VertexInsertion_"+boost::lexical_cast<std::string>(iv+1));
      measurements << alps::ngs::RealVectorObservable("VertexRemoval_"+boost::lexical_cast<std::string>(iv+1));
  }
  measurements << alps::ngs::SimpleRealVectorObservable("MeasurementTimeMsec");
  measurements << alps::ngs::SimpleRealVectorObservable("UpdateTimeMsec");
  measurements << alps::ngs::RealObservable("RecomputeTime");
  for(spin_t flavor=0;flavor<n_flavors;++flavor) {
      std::stringstream tmp;
      tmp<<"VertexHistogram_"<<flavor;
      measurements << alps::ngs::SimpleRealVectorObservable(tmp.str().c_str());
  }
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexInsertion");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexRemoval");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexShift");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexInsertion_count");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexRemoval_count");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexShift_count");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexInsertion_sum");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexRemoval_sum");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexShift_sum");
  measurements << alps::ngs::RealVectorObservable("PerturbationOrderVertex");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexInsertion");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexInsertion_count");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexInsertion_sum");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexRemoval");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexRemoval_count");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexRemoval_sum");
  if (n_multi_vertex_update>1) {
    measurements << alps::ngs::SimpleRealObservable("QuantumNumberConserved");
  }
  measurements << alps::ngs::SimpleRealVectorObservable("PertOrderHistogram");
  if (use_alpha_update)
    measurements << alps::ngs::SimpleRealVectorObservable("AlphaScaleHistogram");
#endif
  measurements.reset(true);
}

///this function is called whenever measurements should be performed. Depending
///on the value of  measurement_method it will choose one particular
///measurement function.
template<class TYPES>
void InteractionExpansion<TYPES>::measure_observables(std::valarray<double>& timings)
{
  assert(timings.size()>=2);
  boost::timer::cpu_timer timer;

  measurements["Sign"] << mycast<REAL_TYPE>(sign);
  if (params.defined("OUTPUT_Sign") ? params["OUTPUT_Sign"] : false) {
      std::cout << " node= " << node << " Sign= " << sign << " pert_order= " << itime_vertices.size() << std::endl;
  }

  pert_order_hist /= pert_order_hist.sum();
  measurements["PertOrderHistogram"] << pert_order_hist;

  const double t1 = timer.elapsed().wall*1E-6;
  if (n_matsubara_measurements>0) {
    compute_W_matsubara();
  }
  const double t2 = timer.elapsed().wall*1E-6;
  if (n_legendre>0) {
    compute_Sl();
  }
  const double t3 = timer.elapsed().wall*1E-6;
  timings[0] = t2-t1;
  timings[1] = t3-t2;

  std::valarray<double> pert_order(n_flavors);
  for(unsigned int i=0;i<n_flavors;++i) { 
    assert(num_rows(M[i].matrix()) == num_cols(M[i].matrix()));
    pert_order[i]=num_rows(M[i].matrix());
  }
  measurements["PertOrder"] << pert_order;

  for(spin_t flavor=0; flavor<n_flavors; ++flavor) {
      std::stringstream tmp;
      tmp<<"VertexHistogram_"<<flavor;
      measurements[tmp.str().c_str()] << vertex_histograms[flavor]->to_valarray();
      vertex_histograms[flavor]->clear();
  }

  if (n_multi_vertex_update>1) {
      measurements["StatisticsVertexInsertion"] << statistics_ins.get_mean();
      measurements["StatisticsVertexRemoval"] << statistics_rem.get_mean();
      measurements["StatisticsDoubleVertexInsertion"] << statistics_dv_ins.get_mean();
      measurements["StatisticsDoubleVertexRemoval"] << statistics_dv_rem.get_mean();

      measurements["StatisticsVertexInsertion_count"] << statistics_ins.get_counter();
      measurements["StatisticsVertexRemoval_count"] << statistics_rem.get_counter();
      measurements["StatisticsDoubleVertexInsertion_count"] << statistics_dv_ins.get_counter();
      measurements["StatisticsDoubleVertexRemoval_count"] << statistics_dv_rem.get_counter();

      measurements["StatisticsVertexInsertion_sum"] << statistics_ins.get_sumval();
      measurements["StatisticsVertexRemoval_sum"] << statistics_rem.get_sumval();
      measurements["StatisticsDoubleVertexInsertion_sum"] << statistics_dv_ins.get_sumval();
      measurements["StatisticsDoubleVertexRemoval_sum"] << statistics_dv_rem.get_sumval();
  }

  measurements["StatisticsVertexShift"] << statistics_shift.get_mean();
  measurements["StatisticsVertexShift_count"] << statistics_shift.get_counter();
  measurements["StatisticsVertexShift_sum"] << statistics_shift.get_sumval();

  statistics_ins.reset();
  statistics_rem.reset();
  statistics_shift.reset();
  statistics_dv_ins.reset();
  statistics_dv_rem.reset();

  for (int iv=0; iv<n_multi_vertex_update; ++iv){
    measurements["VertexInsertion_"+boost::lexical_cast<std::string>(iv+1)] << simple_statistics_ins.get_result(iv);
    measurements["VertexRemoval_"+boost::lexical_cast<std::string>(iv+1)] << simple_statistics_rem.get_result(iv);
  }
  simple_statistics_ins.reset();
  simple_statistics_rem.reset();

  std::valarray<double> pert_vertex(Uijkl.n_vertex_type());
  for (itime_vertex_container::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
    assert(it->type()>=0 && it->type()<Uijkl.n_vertex_type());
    ++pert_vertex[it->type()];
  }
  measurements["PerturbationOrderVertex"] << pert_vertex;
}
