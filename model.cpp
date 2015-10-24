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


double HubbardInteractionExpansion::try_add()
{
  //work array
  std::valarray<int> num_inserted_ops(0,n_flavors);
  //const size_t rank = Uijkl.rank;

  //select a new vertex
  double t = beta*random();
  const size_t v_id = static_cast<size_t>(random()*Uijkl.n_vertex_type());
  const vertex_new new_vertex = Uijkl.get_vertex(v_id);
  const size_t rank = new_vertex.rank();
  assert(rank==2);

  //select a new AF state
  const size_t af_state = static_cast<size_t>(random()*new_vertex.num_af_states());

  const double Uval = new_vertex.Uval();

  std::vector<spin_t> flavor_vertex(rank);
  std::vector<site_t> site_c_vertex(rank), site_a_vertex(rank);
  std::vector<size_t> pos_c_vertex(rank);

  for (size_t i_rank=0; i_rank<rank; ++i_rank) {
    //We insert c^dagger c - alpha.
    size_t flavor_rank = Uijkl.flavor_index(v_id, i_rank);
    size_t site_c = Uijkl.site_index(v_id, 2*i_rank);
    size_t site_a = Uijkl.site_index(v_id, 2*i_rank+1);
    size_t alpha_val = Uijkl.alpha(v_id,af_state,i_rank);

    assert(flavor_rank<n_flavors);
    assert(site_c<n_site);
    assert(site_a<n_site);

    //for later use
    flavor_vertex[i_rank] = flavor_rank;
    site_c_vertex[i_rank] = site_c;
    site_a_vertex[i_rank] = site_a;
    pos_c_vertex[i_rank] = M[flavor_rank].creators().size();

    M[flavor_rank].creators().push_back(creator(flavor_rank, site_c, t, n_matsubara));
    M[flavor_rank].annihilators().push_back(annihilator(flavor_rank, site_a, t, n_matsubara));
    M[flavor_rank].alpha().push_back(alpha_val);

    ++num_inserted_ops[flavor_rank];


  }

  vertices.push_back(vertex(flavor_vertex, site_c_vertex, site_a_vertex, pos_c_vertex, abs_w));
  //vertices.push_back(vertex(flavor0, site, M[0].creators().size()-1, M[0].annihilators().size()-1,
                            //flavor1, site, M[1].creators().size()-1, M[1].annihilators().size()-1, abs_w));

  //int site = (int)(random()*n_site);
  //double alpha0 = random()<0.5?alpha:1-alpha;
  //double alpha1 = 1 - alpha0;
  //spin_t flavor0=0;
  //spin_t flavor1=1;

  //M[0].creators().push_back(creator(flavor0, site, t, n_matsubara));
  //M[0].annihilators().push_back(annihilator(flavor0, site, t, n_matsubara));
  //M[0].alpha().push_back(alpha0); //symmetrized version
//
  //M[1].creators().push_back(creator(flavor1, site, t, n_matsubara));
  //M[1].annihilators().push_back(annihilator(flavor1, site, t,n_matsubara));
  //M[1].alpha().push_back(alpha1); //symmetrized version
  //keep track of vertex list

  //vertices.push_back(vertex(flavor0, site, M[0].creators().size()-1, M[0].annihilators().size()-1,
                        //flavor1, site, M[1].creators().size()-1, M[1].annihilators().size()-1, abs_w));

  //perform fastupdate up for weight
  double lambda0=fastupdate_up(flavor0, true); // true means compute_only_weight
  double lambda1=fastupdate_up(flavor1, true);
  double num_vertex_types=n_site;
  double metropolis_weight=-beta*abs_w*num_vertex_types/(vertices.size())*lambda0*lambda1;
  //return weight
  return metropolis_weight;
  return 0.0;
}


void HubbardInteractionExpansion::perform_add()
{
  //perform the fastupdate up move
  fastupdate_up(0,false);
  fastupdate_up(1,false);
}


void HubbardInteractionExpansion::reject_add()
{
  //get rid of the operators
  M[0].creators().pop_back();
  M[0].annihilators().pop_back();
  M[0].alpha().pop_back();
  M[1].creators().pop_back();
  M[1].annihilators().pop_back();
  M[1].alpha().pop_back();
  //get rid of the vertex from vertex list
  vertices.pop_back();
}


double HubbardInteractionExpansion::try_remove(unsigned int vertex_nr)
{
  //get weight
  double lambda0 = fastupdate_down(vertex_nr, 0, true);  // true means compute_only_weight
  double lambda1 = fastupdate_down(vertex_nr, 1, true);  
  double pert_order=num_rows(M[0].matrix());
  double num_vertex_types=n_site;
  //return weight
  return  -pert_order/(beta*onsite_U*num_vertex_types)*lambda0*lambda1;
}


void HubbardInteractionExpansion::perform_remove(unsigned int vertex_nr)
{
  //perform fastupdate down
  fastupdate_down(vertex_nr, 0, false);  // false means really perform, not only compute weight
  fastupdate_down(vertex_nr, 1, false);  // false means really perform, not only compute weight
  //get rid of operators
  M[0].creators().pop_back();
  M[0].annihilators().pop_back();
  M[0].alpha().pop_back();
  M[1].creators().pop_back();
  M[1].annihilators().pop_back();
  M[1].alpha().pop_back();
  //get rid of vertex list entries
  vertices.pop_back();
}


void HubbardInteractionExpansion::reject_remove()
{
  //do nothing
  return;
}
