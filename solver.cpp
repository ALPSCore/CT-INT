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

#include <alps/numeric/matrix/algorithms.hpp>

///The basic updates for the InteractionExpansion algorithm: adding and removing vertices.
///This is the heart of InteractionExpansion's code.
void InteractionExpansion::interaction_expansion_step(void)
{
  //temp work memory
  static fastupdate_add_helper add_helper(n_flavors);
  static fastupdate_remove_helper remove_helper(n_flavors);

  size_t nv_updated = (size_t) (random()*n_multi_vertex_update)+1;

//#ifndef NDEBUG
  /***** VERY EXPENSIVE TEST *****/
  //double det_updated = 1/M.determinant();
//#endif

  const int pert_order= itime_vertices.size();   //current order of perturbation series
  double metropolis_weight=0.;
  double det_rat=0;
  //static unsigned int i=0; ++i;
  if(random()<0.5){  //trying to ADD vertex
    M.sanity_check(itime_vertices);
    if(itime_vertices.size()>=max_order)
      return; //we have already reached the highest perturbation order
    std::vector<itime_vertex> new_vertices = generate_vertices(Uijkl,random,beta,nv_updated);
    if (!is_quantum_number_conserved(new_vertices)) {
      simple_statistics_ins.not_valid_state(nv_updated-1);
      return;
    }
    if (!is_irreducible(new_vertices)) {
      simple_statistics_ins.reducible(nv_updated-1);
      return;
    }

    boost::tie(metropolis_weight,det_rat)=try_add(add_helper,nv_updated, new_vertices);
    if (nv_updated>=2) {
      statistics_ins.add_sample(compute_spread(new_vertices,beta), std::min(fabs(metropolis_weight),1.0), nv_updated-2);
    }
    if(fabs(metropolis_weight)> random()){
      //measurements["VertexInsertion"]<<1.;
      perform_add(add_helper,nv_updated);
      sign*=metropolis_weight<0?-1:1;
      M.sanity_check(itime_vertices);
      assert(is_quantum_number_conserved(new_vertices));
      simple_statistics_ins.accepted(nv_updated-1);
    }else{
      //measurements["VertexInsertion"]<<0.;
      reject_add(add_helper,nv_updated);
      M.sanity_check(itime_vertices);
      simple_statistics_ins.rejected(nv_updated-1);
    }
  }else{ // try to REMOVE a vertex
    M.sanity_check(itime_vertices);
    if(pert_order < nv_updated) {
      return;
    }

    //choose vertices to be removed
    const std::vector<size_t> vertices_nr = pickup_a_few_numbers(pert_order, nv_updated, random);
    std::vector<itime_vertex> vertices_to_be_removed(nv_updated);
    for (int iv=0; iv<nv_updated; ++iv) {
      vertices_to_be_removed[iv] = itime_vertices[vertices_nr[iv]];
    }
    if (!is_quantum_number_conserved(vertices_to_be_removed)) {
      simple_statistics_rem.not_valid_state(nv_updated-1);
      return;
    }
    if (!is_irreducible(vertices_to_be_removed)) {
      simple_statistics_rem.reducible(nv_updated - 1);
      return;
    }

    boost::tie(metropolis_weight,det_rat)=try_remove(vertices_nr, remove_helper); //get the determinant ratio. don't perform fastupdate yet
    if(fabs(metropolis_weight)> random()){ //do the actual update
      //measurements["VertexRemoval"]<<1.;
      perform_remove(vertices_nr, remove_helper);
      sign*=metropolis_weight<0?-1:1;
      M.sanity_check(itime_vertices);
      simple_statistics_rem.accepted(nv_updated-1);
    }else{
      //measurements["VertexRemoval"]<<0.;
      reject_remove(remove_helper);
      M.sanity_check(itime_vertices);
      simple_statistics_rem.rejected(nv_updated-1);
    }
  }//end REMOVE
  weight=metropolis_weight;
  for(spin_t flavor=0; flavor<n_flavors; ++flavor) {
    vertex_histograms[flavor]->count(num_rows(M[flavor].matrix()));
  }
}

///Every now and then we have to recreate M from scratch to avoid roundoff
///error. This is done by iserting the vertices starting from zero.
void InteractionExpansion::reset_perturbation_series()
{
  //static fastupdate_add_helper add_helper(n_flavors);
  big_inverse_m_matrix M2(M); //make a copy of M

  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    alps::numeric::matrix<GTYPE> G0(M[flavor].matrix());
    std::fill(G0.get_values().begin(), G0.get_values().end(), 0);

    const size_t Nv = M[flavor].matrix().num_rows();

    if (Nv==0) {
      M[flavor].matrix() = alps::numeric::matrix<GTYPE>(0,0);
    } else {
      //construct G0 matrix
      for (size_t q=0; q<Nv; ++q) {
        for (size_t p=0; p<Nv; ++p) {
          G0(p, q) = green0_spline_for_M(flavor, p, q);
        }
      }
      for (size_t p=0; p<Nv; ++p) {
        G0(p, p) -= M[flavor].alpha()[p];
      }
      M[flavor].matrix() = alps::numeric::inverse(G0);
    }
  }

  for(unsigned int z=0;z<M2.size();++z){
    double max_diff=0;
    for(unsigned int j=0;j<num_cols(M2[z].matrix());++j){
      for(unsigned int i=0;i<num_rows(M2[z].matrix());++i){
        double diff=M[z].matrix()(i,j)-M2[z].matrix()(i,j);
        if(std::abs(diff)>max_diff) max_diff=std::abs(diff);
      }
    }
    //std::cout << "debug: max_diff = " << max_diff << std::endl;
    if(max_diff > 1.e-8)
      std::cout<<"WARNING: roundoff errors in flavor: "<<z<<" max diff "<<max_diff<<std::endl;
  }

#ifndef NDEBUG
  //recompute weight
#endif
}

