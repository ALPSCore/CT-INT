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

  int pert_order=vertices_new.size();   //current order of perturbation series
  double metropolis_weight=0.;
  static unsigned int i=0; ++i;
  if(random()<0.5){  //trying to ADD vertex
    if(vertices_new.size()>=max_order)
      return; //we have already reached the highest perturbation order
    metropolis_weight=try_add(add_helper,1);
    if(fabs(metropolis_weight)> random()){
      measurements["VertexInsertion"]<<1.;
      perform_add(add_helper,1);
      sign*=metropolis_weight<0?-1:1;
      assert(num_rows(M[0].matrix())==vertices_new.size());
    }else{
      measurements["VertexInsertion"]<<0.;
      reject_add(add_helper,1);
      assert(num_rows(M[0].matrix())==vertices_new.size());
    }
  }else{ // try to REMOVE a vertex
    pert_order=vertices_new.size(); //choose a vertex
    if(pert_order < 1) {
      return;     //we have an empty list or one with just one vertex
    }
    //this might be the ideal place to do some cleanup, e.g. get rid of the roundoff errors and such.
    int vertex_nr=(int)(random() * pert_order);
    metropolis_weight=try_remove(vertex_nr, remove_helper); //get the determinant ratio. don't perform fastupdate yet
    if(fabs(metropolis_weight)> random()){ //do the actual update
      measurements["VertexRemoval"]<<1.;
      perform_remove(vertex_nr, remove_helper);
      sign*=metropolis_weight<0?-1:1;
    }else{
      measurements["VertexRemoval"]<<0.;
      reject_remove(remove_helper);
    }
  }//end REMOVE

  /**** sanity check for onsite U ******/
  M.sanity_check();
  {
    size_t Nv = 0;
    for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
      Nv += M[flavor].creators().size();
    }
    assert(2*vertices_new.size()==Nv);
  }

  weight=metropolis_weight;
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

    //construct G0 matrix
    for (size_t q=0; q<Nv; ++q) {
      for (size_t p=0; p<Nv; ++p) {
        G0(p, q) = green0_spline_new(M[flavor].annihilators()[p], M[flavor].creators()[q]);
      }
    }
    for (size_t p=0; p<Nv; ++p) {
      G0(p, p) += M[flavor].alpha()[p];
    }

    M[flavor].matrix() = alps::numeric::inverse(G0);
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
}

