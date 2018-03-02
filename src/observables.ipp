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

#include "accumulators.hpp"
#include <complex>

///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables. It is also called at the start of every DMFT
///iteration.
template<class TYPES>
void InteractionExpansion<TYPES>::initialize_observables(void)
{
  if(measurements.has("Sign")){
    measurements.clear();
  }
  measurements << SimpleRealObservable("Sign");
  measurements << SimpleRealVectorObservable("PertOrder");
  if(n_matsubara_measurements>0) {
    for(unsigned int flavor=0;flavor<n_flavors;++flavor){
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

  if (n_legendre>0) {
    for(unsigned int flavor=0;flavor<n_flavors;++flavor){
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
  for(int flavor=0;flavor<n_flavors;++flavor) {
    measurements << SimpleRealVectorObservable("densities_"+boost::lexical_cast<std::string>(flavor));
  }
  measurements << SimpleRealObservable("density_correlation");
  measurements << SimpleRealVectorObservable("n_i n_j");

  for(unsigned int flavor=0;flavor<n_flavors;++flavor){
    for(unsigned int i=0;i<n_site;++i){
      std::stringstream density_name, sz_name;
      density_name<<"density_"<<flavor;
      if (n_site>1) density_name<<"_"<<i;
      measurements << SimpleRealObservable(density_name.str().c_str());
    }
  }
  for(unsigned int i=0;i<n_site;++i){
    std::stringstream sz_name, sz2_name, sz0_szj_name;
    sz_name<<"Sz_"<<i;
    sz2_name<<"Sz2_"<<i;
    sz0_szj_name<<"Sz0_Sz"<<i;
  }
  /*
  for (int iv=0; iv<n_multi_vertex_update; ++iv)  {
      measurements << SimpleRealVectorObservable("VertexInsertion_"+boost::lexical_cast<std::string>(iv+1));
      measurements << SimpleRealVectorObservable("VertexRemoval_"+boost::lexical_cast<std::string>(iv+1));
  }
  */
  measurements << SimpleRealVectorObservable("MeasurementTimeMsec");
  measurements << SimpleRealVectorObservable("UpdateTimeMsec");
  measurements << SimpleRealVectorObservable("UpdateTimeMsecAllWalkers");
  measurements << SimpleRealObservable("RecomputeTime");
  for(spin_t flavor=0;flavor<n_flavors;++flavor) {
      std::stringstream tmp;
      tmp<<"VertexHistogram_"<<flavor;
      measurements << SimpleRealVectorObservable(tmp.str().c_str());
  }

  measurements << SimpleRealVectorObservable("PerturbationOrderVertex");
  measurements << SimpleRealVectorObservable("PertOrderHistogram");
  measurements << SimpleRealVectorObservable("ACCEPTANCE_RATE_EXCHANGE");
  //measurements.reset(true);
}

///this function is called whenever measurements should be performed. Depending
///on the value of  measurement_method it will choose one particular
///measurement function.
template<class TYPES>
void InteractionExpansion<TYPES>::measure_observables(std::valarray<double>& timings)
{
  assert(timings.size()>=2);
  boost::timer::cpu_timer timer;

  //compute M from A
  submatrix_update->compute_M(M_flavors);
  const M_TYPE sign = submatrix_update->sign();

  measurements["Sign"] << mycast<REAL_TYPE>(sign);
  if (parms.defined("OUTPUT_Sign") ? parms["OUTPUT_Sign"] : false) {
      std::cout << " node= " << node << " Sign= " << sign << " pert_order= " << submatrix_update->pert_order() << std::endl;
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
    pert_order[i]=M_flavors[i].size1();
  }
  measurements["PertOrder"] << pert_order;

  for(spin_t flavor=0; flavor<n_flavors; ++flavor) {
      std::stringstream tmp;
      tmp<<"VertexHistogram_"<<flavor;
      measurements[tmp.str().c_str()] << vertex_histograms[flavor]->to_valarray();
      vertex_histograms[flavor]->clear();
  }

  std::valarray<double> pert_vertex(Uijkl.n_vertex_type());
  const itime_vertex_container& itime_vertices = submatrix_update->itime_vertices();
  for (itime_vertex_container::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
    assert(it->type()>=0 && it->type()<Uijkl.n_vertex_type());
    ++pert_vertex[it->type()];
  }
  measurements["PerturbationOrderVertex"] << pert_vertex;
}
