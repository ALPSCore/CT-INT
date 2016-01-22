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

///The basic updates for the InteractionExpansion algorithm: adding and removing vertices.
///This is the heart of InteractionExpansion's code.
template<class TYPES>
void InteractionExpansion<TYPES>::removal_insertion_update(void) {
  const int nv_updated = update_prop.gen_Nv(random.engine());

  //boost::timer::cpu_timer timer;
  if (nv_updated==1) {
    removal_insertion_single_vertex_update();
  } else if (nv_updated==2) {
    removal_insertion_double_vertex_update();
  } else {
    multi_vertex_update(nv_updated);
  }

  for(spin_t flavor=0; flavor<n_flavors; ++flavor) {
    vertex_histograms[flavor]->count(num_rows(M[flavor].matrix()));
  }
}

template<class TYPES>
void InteractionExpansion<TYPES>::removal_insertion_single_vertex_update(void)
{
  const int nv_updated = 1;

  const int pert_order= itime_vertices.size();   //current order of perturbation series
  M_TYPE metropolis_weight=0.;
  M_TYPE det_rat=0;
  if(random()<0.5){  //trying to ADD vertex
    boost::timer::cpu_timer timer;
    M.sanity_check(itime_vertices);
    if(pert_order+nv_updated>max_order)
      return; //we have already reached the highest perturbation order

    std::vector<itime_vertex> new_vertices;
    std::pair<int,int> v_pair;
    if (single_vertex_update_non_density_type) {
      new_vertices = generate_itime_vertices(Uijkl,random,beta,nv_updated,all_type());
    } else {
      new_vertices = generate_itime_vertices(Uijkl,random,beta,nv_updated,density_type());
    }

    assert(new_vertices.size()==nv_updated || new_vertices.size()==0);
    if (new_vertices.size()<nv_updated) {
      simple_statistics_ins.not_valid_state(nv_updated-1);
      update_prop.generated_invalid_update(nv_updated);
      return;
    }

    itime_vertex_container itime_vertices_new(itime_vertices);
    for (itime_vertex_container::const_iterator it=new_vertices.begin(); it!=new_vertices.end(); ++it) {
      itime_vertices_new.push_back(*it);
    }
    if (force_quantum_number_conservation &&
        (!is_quantum_number_conserved(itime_vertices_new) || !is_quantum_number_within_range(itime_vertices_new))) {
      return;
    }

    update_prop.generated_valid_update(nv_updated);
    boost::tie(metropolis_weight,det_rat)=try_add(nv_updated, new_vertices);

    double p_ins, p_rem;
    if (single_vertex_update_non_density_type) {
      p_ins = 1.0 / (beta * Uijkl.n_vertex_type());
      p_rem = 1.0 / itime_vertices_new.size();
    } else {
      p_ins = 1.0 / (beta * Uijkl.num_density_vertex_type());
      p_rem = 1.0 / std::count_if(itime_vertices_new.begin(),itime_vertices_new.end(),density_type());
    }
    metropolis_weight *= p_rem/p_ins;

    if(std::abs(metropolis_weight)> random()){
      perform_add(nv_updated);
      sign*=metropolis_weight/std::abs(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      simple_statistics_ins.accepted(nv_updated-1);
      if(std::abs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);
      //std::cout << "add accepted " << timer.elapsed().wall*1E-6 << std::endl;
    }else{
      reject_add(nv_updated);
      M.sanity_check(itime_vertices);
      simple_statistics_ins.rejected(nv_updated-1);
      //std::cout << "add rejected " << timer.elapsed().wall*1E-6 << std::endl;
    }
#ifndef NDEBUG
    sanity_check();
#endif
  }else{ // try to REMOVE a vertex
    M.sanity_check(itime_vertices);
    if(pert_order < nv_updated) {
      return;
    }

    //choose vertices to be removed
    std::vector<int> vertices_nr;
    int n_vpair;
    double F;
    double p_ins, p_rem;
    if (single_vertex_update_non_density_type) {
      vertices_nr = pick_up_itime_vertices(itime_vertices, random, nv_updated, all_type());
    } else {
      vertices_nr = pick_up_itime_vertices(itime_vertices, random, nv_updated, density_type());
    }

    if (vertices_nr.size()==0)
      return;

    itime_vertex_container vertices_to_be_removed(nv_updated), itime_vertices_new(itime_vertices);
    for (int iv=0; iv<nv_updated; ++iv) {
      assert(vertices_nr[iv]<itime_vertices.size());
      vertices_to_be_removed[iv] = itime_vertices[vertices_nr[iv]];
    }
    remove_elements_from_vector(itime_vertices_new, vertices_nr);

    if (force_quantum_number_conservation &&
        (!is_quantum_number_conserved(itime_vertices_new) || !is_quantum_number_within_range(itime_vertices_new))) {
      return;
    }

    if (single_vertex_update_non_density_type) {
      p_ins = 1.0/(beta*Uijkl.n_vertex_type());
      p_rem = 1.0/pert_order;
    } else {
      p_ins = 1.0/(beta*Uijkl.num_density_vertex_type());
      p_rem = 1.0/std::count_if(itime_vertices.begin(), itime_vertices.end(), density_type());
    }

    update_prop.generated_valid_update(nv_updated);

    boost::tie(metropolis_weight,det_rat)=try_remove(vertices_nr); //get the determinant ratio. don't perform fastupdate yet
    metropolis_weight *= p_ins/p_rem;
    if(std::abs(metropolis_weight)> random()){ //do the actual update
      perform_remove(vertices_nr);
      sign*=metropolis_weight/std::abs(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      simple_statistics_rem.accepted(nv_updated-1);
      if(std::abs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);
    }else{
      reject_remove();
      M.sanity_check(itime_vertices);
      simple_statistics_rem.rejected(nv_updated-1);
    }
#ifndef NDEBUG
    sanity_check();
#endif
  }//end REMOVE
  weight=std::abs(metropolis_weight);
}

template<class TYPES>
void InteractionExpansion<TYPES>::removal_insertion_double_vertex_update(void)
{
  const int nv_updated = 2;

  if (mv_update_valid_pair.size()==0)
    return;

  const int pert_order= itime_vertices.size();   //current order of perturbation series
  M_TYPE metropolis_weight=0.;
  M_TYPE det_rat=0;
  if(random()<0.5){  //trying to ADD vertex
    M.sanity_check(itime_vertices);
    if(pert_order+nv_updated>max_order)
      return; //we have already reached the highest perturbation order

    add_helper.op = non_density_type_in_window(beta*random(),
                                               std::min(beta,window_width+window_dist(random.engine())),
                                               beta);

    std::vector<itime_vertex> new_vertices;
    std::pair<int,int> v_pair;
    v_pair = mv_update_valid_pair[mv_update_valid_pair.size()*random()];
    new_vertices = generate_valid_vertex_pair2(Uijkl,v_pair,random,beta,symm_exp_dist);

    assert(new_vertices.size()==nv_updated || new_vertices.size()==0);
    if (new_vertices.size()<nv_updated) {
      simple_statistics_ins.not_valid_state(nv_updated-1);
      update_prop.generated_invalid_update(nv_updated);
      return;
    }

    itime_vertex_container itime_vertices_new(itime_vertices);
    for (itime_vertex_container::const_iterator it=new_vertices.begin(); it!=new_vertices.end(); ++it) {
      itime_vertices_new.push_back(*it);
    }
    if (force_quantum_number_conservation &&
        (!is_quantum_number_conserved(itime_vertices_new) || !is_quantum_number_within_range(itime_vertices_new))) {
      return;
    }

    update_prop.generated_valid_update(nv_updated);
    boost::tie(metropolis_weight,det_rat)=try_add(nv_updated, new_vertices);

    //double p_ins, p_rem;
    const double dtau = mymod(new_vertices[1].time()-new_vertices[0].time(), beta);
    int n_vpair;
    double F;
    pick_up_valid_vertex_pair2(itime_vertices_new,
                                 std::make_pair(new_vertices[0].type(),new_vertices[1].type()),
                                 beta, symm_exp_dist, random, n_vpair, F);
    if (n_vpair==0)
      throw std::logic_error("v_pair must be larger than 0.");
    metropolis_weight *= (beta*beta)*symm_exp_dist.coeff_X(dtau)/F;

    statistics_ins.add_sample(compute_spread(new_vertices,beta), std::min(std::abs(metropolis_weight),1.0), nv_updated-2);
    statistics_dv_ins.add_sample(mymod(new_vertices[1].time()-new_vertices[0].time(), beta),
                                 std::min(std::abs(metropolis_weight), 1.0),
                                 std::distance(mv_update_valid_pair.begin(), std::find(mv_update_valid_pair.begin(), mv_update_valid_pair.end(), v_pair))
    );
    if(std::abs(metropolis_weight)> random()){
      perform_add(nv_updated);
      sign*=metropolis_weight/std::abs(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      simple_statistics_ins.accepted(nv_updated-1);
      if(std::abs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);
    }else{
      reject_add(nv_updated);
      M.sanity_check(itime_vertices);
      simple_statistics_ins.rejected(nv_updated-1);
    }
#ifndef NDEBUG
    sanity_check();
#endif
  }else{ // try to REMOVE a vertex
    M.sanity_check(itime_vertices);
    if(pert_order < nv_updated) {
      return;
    }
    remove_helper.op = non_density_type_in_window(beta*random(),
                                                  std::min(beta,window_width+window_dist(random.engine())),
                                                  beta);

    //choose vertices to be removed
    std::vector<int> vertices_nr;
    int n_vpair;
    double F;
    std::pair<int,int> v_pair = mv_update_valid_pair[mv_update_valid_pair.size()*random()];
    std::pair<int,int> r = pick_up_valid_vertex_pair2(itime_vertices, v_pair, beta, symm_exp_dist, random, n_vpair, F);
    if (n_vpair>0) {
      vertices_nr.resize(2);
      vertices_nr[0]=r.first;
      vertices_nr[1]=r.second;
    } else {
      vertices_nr.resize(0);
    }

    if (vertices_nr.size()==0)
      return;

    std::vector<itime_vertex> vertices_to_be_removed(nv_updated);
    itime_vertex_container itime_vertices_new(itime_vertices);
    for (int iv=0; iv<nv_updated; ++iv) {
      assert(vertices_nr[iv]<itime_vertices.size());
      vertices_to_be_removed[iv] = itime_vertices[vertices_nr[iv]];
    }
    remove_elements_from_vector(itime_vertices_new, vertices_nr);

    if (force_quantum_number_conservation &&
        (!is_quantum_number_conserved(itime_vertices_new) || !is_quantum_number_within_range(itime_vertices_new))) {
      return;
    }

    assert(vertices_to_be_removed[1].type()>vertices_to_be_removed[0].type());

    update_prop.generated_valid_update(nv_updated);

    boost::tie(metropolis_weight,det_rat)=try_remove(vertices_nr); //get the determinant ratio. don't perform fastupdate yet
    metropolis_weight /= (beta*beta)*
                         symm_exp_dist.coeff_X(
                             mymod(vertices_to_be_removed[1].time()-vertices_to_be_removed[0].time(), beta)
                         )/F;
    statistics_rem.add_sample(compute_spread(vertices_to_be_removed, beta), std::min(std::abs(metropolis_weight), 1.0),
                              nv_updated - 2);
    int idx = std::distance(mv_update_valid_pair.begin(),
                            std::find(
                                mv_update_valid_pair.begin(),
                                mv_update_valid_pair.end(),
                                std::make_pair(vertices_to_be_removed[0].type(),vertices_to_be_removed[1].type())
                            )
    );
    statistics_dv_rem.add_sample(
        mymod(vertices_to_be_removed[1].time()-vertices_to_be_removed[0].time(), beta),
        std::min(std::abs(metropolis_weight), 1.0),
        idx);
    if(std::abs(metropolis_weight)> random()){ //do the actual update
      perform_remove(vertices_nr);
      sign*=metropolis_weight/std::abs(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      simple_statistics_rem.accepted(nv_updated-1);
      if(std::abs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);
    }else{
      reject_remove();
      M.sanity_check(itime_vertices);
      simple_statistics_rem.rejected(nv_updated-1);
    }
#ifndef NDEBUG
    sanity_check();
#endif
  }//end REMOVE
  weight=std::abs(metropolis_weight);
}

template<class TYPES>
void InteractionExpansion<TYPES>::multi_vertex_update(int nv_updated)
{
  const int pert_order= itime_vertices.size();   //current order of perturbation series
  M_TYPE metropolis_weight=0.;
  M_TYPE det_rat=0;

  if (nv_updated<2) {
    throw std::logic_error("Error: nv_updated<2.");
  }

  if(random()<0.5){  //trying to ADD vertex
    M.sanity_check(itime_vertices);
    if(pert_order+nv_updated>max_order) {
      return; //we have already reached the highest perturbation order
    }

    add_helper.op = non_density_type_in_window(beta*random(),
                                               std::min(beta,window_width+window_dist(random.engine())),
                                               beta);
    std::vector<itime_vertex> new_vertices =
      generate_itime_vertices(Uijkl,random,beta,nv_updated,add_helper.op);

    assert(new_vertices.size()==nv_updated || new_vertices.size()==0);
    itime_vertex_container itime_vertices_new(itime_vertices);
    for (std::vector<itime_vertex>::const_iterator it=new_vertices.begin(); it!=new_vertices.end(); ++it) {
      itime_vertices_new.push_back(*it);
    }
    if (!is_quantum_number_conserved(itime_vertices_new)) {
      simple_statistics_ins.not_valid_state(nv_updated-1);
      update_prop.generated_invalid_update(nv_updated);
      return;
    }

    update_prop.generated_valid_update(nv_updated);

    boost::tie(metropolis_weight,det_rat)=try_add(nv_updated, new_vertices);
    statistics_ins.add_sample(compute_spread(new_vertices,beta), std::min(std::abs(metropolis_weight),1.0), nv_updated-2);
    //std::cout << "mw " << metropolis_weight << std::endl;

    if(std::abs(metropolis_weight)> random()){
      //std::cout << "accepted ins " << std::endl;
      perform_add(nv_updated);
      //std::cout << "accepted ins (2)" << std::endl;
      sign*=metropolis_weight/std::abs(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      //std::cout << "accepted ins (3)" << std::endl;
      //std::cout << "debug1 " << simple_statistics_ins.counter_[nv_updated-1][2] << std::endl;
      //std::cout << "debug1 " << simple_statistics_ins.counter_[nv_updated-1][3] << std::endl;
      simple_statistics_ins.accepted(nv_updated-1);
      //std::cout << "debug2 " << simple_statistics_ins.counter_[nv_updated-1][2] << std::endl;
      //std::cout << "debug2 " << simple_statistics_ins.counter_[nv_updated-1][3] << std::endl;
      if(std::abs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);
    }else{
      //std::cout << "rejected ins " << std::endl;
      reject_add(nv_updated);
      M.sanity_check(itime_vertices);
      simple_statistics_ins.rejected(nv_updated-1);
    }
#ifndef NDEBUG
    sanity_check();
#endif
  }else{ // try to REMOVE a vertex
    M.sanity_check(itime_vertices);
    if(pert_order < nv_updated) {
      return;
    }
    remove_helper.op = non_density_type_in_window(beta*random(),
            std::min(beta,window_width+window_dist(random.engine())),
            beta);

    //choose vertices to be removed
    const std::vector<int>& vertices_nr = pick_up_itime_vertices(itime_vertices, random, nv_updated, remove_helper.op);
    if (vertices_nr.size()==0) {
      return;
    }
    itime_vertex_container vertices_to_be_removed(nv_updated), itime_vertices_new(itime_vertices);
    for (int iv=0; iv<nv_updated; ++iv) {
      vertices_to_be_removed[iv] = itime_vertices[vertices_nr[iv]];
    }
    remove_elements_from_vector(itime_vertices_new, vertices_nr);
    //std::cout << "debug " << is_quantum_number_conserved(itime_vertices_new) << std::endl;
    if (!is_quantum_number_conserved(itime_vertices_new)) {
      simple_statistics_rem.not_valid_state(nv_updated-1);
      update_prop.generated_invalid_update(nv_updated);
      return;
    }

    update_prop.generated_valid_update(nv_updated);

    boost::tie(metropolis_weight,det_rat)=try_remove(vertices_nr); //get the determinant ratio. don't perform fastupdate yet
    //std::cout << " mw " << metropolis_weight << std::endl;
    statistics_rem.add_sample(compute_spread(vertices_to_be_removed, beta), std::min(std::abs(metropolis_weight), 1.0),
                              nv_updated - 2);
    if(std::abs(metropolis_weight)> random()){ //do the actual update
      //std::cout << "accepted rem " << std::endl;
      perform_remove(vertices_nr);
      sign*=metropolis_weight/std::abs(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      simple_statistics_rem.accepted(nv_updated-1);

      if(std::abs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);
    }else{
      //std::cout << "rejected rem " << std::endl;
      reject_remove();
      M.sanity_check(itime_vertices);
      simple_statistics_rem.rejected(nv_updated-1);
    }
#ifndef NDEBUG
    sanity_check();
#endif
  }//end REMOVE
  weight=metropolis_weight;
}

//Shift updates: shift a vertex.
template<class TYPES>
void InteractionExpansion<TYPES>::shift_update(void) {
  const int pert_order= itime_vertices.size();//CHECKED
  if (pert_order<=1)
    return;

  //choose a truely-non-density-type vertex
  std::vector<int> vertices_nr;
  vertices_nr.reserve(pert_order);
  for (int iv=0; iv<itime_vertices.size(); ++iv) {
    if (Uijkl.get_is_truely_non_density_type()[itime_vertices[iv].type()])
      vertices_nr.push_back(iv);
  }
  if (vertices_nr.size()==0)
    return;

  const int iv = vertices_nr[vertices_nr.size()*random()];
  //std::cout << " debug shift " << itime_vertices[iv].type() << std::endl;
  const double new_time = shift_helper.new_itime(itime_vertices[iv].time(), beta, random.engine());
  const double diff_time = std::abs(new_time-itime_vertices[iv].time());

  //check quantum number
  if (force_quantum_number_conservation) {
    const double old_time = itime_vertices[iv].time();
    itime_vertices[iv].set_time(new_time);

    if (is_quantum_number_conserved(itime_vertices) && is_quantum_number_within_range(itime_vertices)) {
      itime_vertices[iv].set_time(old_time);
    } else {
      itime_vertices[iv].set_time(old_time);
      //std::cout << "rejected" << std::endl;
      return;
    }
  }

  M_TYPE metropolis_weight = try_shift(iv, new_time);

  statistics_shift.add_sample(std::min(diff_time,beta-diff_time), std::min(std::abs(metropolis_weight),1.0), 0);
  if(std::abs(metropolis_weight)> random()){ //do the actual update
    perform_shift(iv);
    sign*=metropolis_weight/std::abs(metropolis_weight);
    det*=metropolis_weight;
    if(std::abs(metropolis_weight)<1E-5)
      reset_perturbation_series(false);
#ifndef NDEBUG
    M.sanity_check(itime_vertices);
#endif
    //std::cout << "Done shift " << metropolis_weight << std::endl;
  }else {
    reject_shift(iv);
#ifndef NDEBUG
    M.sanity_check(itime_vertices);
#endif
  }
#ifndef NDEBUG
  sanity_check();
#endif
}

//Spin flip updates: change the spin of a vertex
template<class TYPES>
void InteractionExpansion<TYPES>::spin_flip_update(void) {
  const int pert_order = itime_vertices.size();
  if (pert_order==0) {
    return;
  }

  //choose a vertex and a new af state
  const int iv = (int) pert_order*random();
  itime_vertex& v = itime_vertices[iv];
  const vertex_definition<M_TYPE> v_def = Uijkl.get_vertex(v.type());
  if (Uijkl.get_vertex(v.type()).num_af_states()==1) {
    return;
  }
  int new_af_state;
  const int old_af_state = v.af_state();
  while (true) {
    new_af_state = (int) Uijkl.get_vertex(v.type()).num_af_states()*random();
    if (new_af_state!=v.af_state()) {
      break;
    }
  }

  //std::cout << "debug v.type " << v.type() << std::endl;
  //std::cout << "old_af " << old_af_state << std::endl;
  //std::cout << "new_af " << new_af_state << std::endl;
  //std::cout << "old_M " << M[0].matrix() << std::endl;

  //change status
  v.set_af_state(new_af_state);
  std::vector<std::vector<int> > positions(n_flavors);
  std::vector<std::vector<M_TYPE> > old_alpha(n_flavors);
  for (int i_rank=0; i_rank<v_def.rank(); ++i_rank) {
    const int flavor_rank = v_def.flavors()[i_rank];
    int pos = M[flavor_rank].find_row_col(v.time(), v.type(), i_rank);
    positions[flavor_rank].push_back(pos);
    old_alpha[flavor_rank].push_back(M[flavor_rank].bare_alpha_at(pos));
    //std::cout << "debug i_rank " << i_rank << " flavor " << flavor_rank << " " << M[flavor_rank].bare_alpha_at(pos) << " new " << v_def.get_alpha(new_af_state,i_rank) << std::endl;
    if (M[flavor_rank].bare_alpha_at(pos)!=v_def.get_alpha(new_af_state,i_rank)) {
      M[flavor_rank].set_alpha(pos, v_def.get_alpha(new_af_state,i_rank));
    }
  }

  //update M
  M_TYPE det_rat = 1.0;
  for (int flavor=0; flavor<n_flavors; ++flavor) {
    if (positions[flavor].size()>0) {
      //std::cout << "debug flavor " << flavor << std::endl;
      assert(positions[flavor].size()==old_alpha[flavor].size());
      det_rat *= fastupdate_spin_flip(flavor, positions[flavor], old_alpha[flavor], true);
    }
  }

  //std::cout << "debug, det_rat " << det_rat << std::endl;

  if(std::abs(det_rat)> random()) { //do the actual update
    //std::cout << "accepted " << std::endl;
    sign *= det_rat/std::abs(det_rat);
    det *= det_rat;
    for (int flavor=0; flavor<n_flavors; ++flavor) {
      if (positions[flavor].size()>0) {
        fastupdate_spin_flip(flavor, positions[flavor], old_alpha[flavor], false);
      }
    }
    //std::cout << "new_M " << M[0].matrix() << std::endl;
  } else {
    //std::cout << "rejected " << std::endl;
    itime_vertices[iv].set_af_state(old_af_state);
    for (int flavor=0; flavor<n_flavors; ++flavor) {
      for (int i_op=0; i_op<positions[flavor].size(); ++i_op) {
        M[flavor].set_alpha(positions[flavor][i_op], old_alpha[flavor][i_op]);
      }
    }
  }
  M.sanity_check(itime_vertices);
  sanity_check();
}

//Update alpha for non-density-type vertices
template<class TYPES>
void InteractionExpansion<TYPES>::alpha_update(void) {
  if (itime_vertices.size()==0 || alpha_scale_max==1.0)
    return;

  const int new_alpha_scale_idx = static_cast<int>(num_alpha_scale_values*random());
  const double new_alpha_scale = alpha_scale_values[new_alpha_scale_idx];

  const M_TYPE det_old = det;
  const M_TYPE sign_old = sign;
  const big_inverse_m_matrix<M_TYPE> M_old(M);

  M.set_alpha_scale(new_alpha_scale);
  boost::timer::cpu_timer timer;
  reset_perturbation_series(false);
  const M_TYPE metropolis_weight = (det/det_old) * flat_histogram_alpha.weight_ratio(new_alpha_scale_idx, alpha_scale_idx);

  if(std::abs(metropolis_weight)> random()){ //do the actual update
    alpha_scale_idx = new_alpha_scale_idx;
  }else {
    det = det_old;
    M = M_old;
    sign = sign_old;
  }
}

///We recreate M from scratch to avoid roundoff error.
//This is done by constructing a matrix consisting of bare Green's function and invert it.
template<class TYPES>
void InteractionExpansion<TYPES>::reset_perturbation_series(bool verbose)
{
  typedef alps::numeric::matrix<M_TYPE> matrix_t;
  big_inverse_m_matrix<M_TYPE> M2;
  M_TYPE det_old;
  if (verbose) {
    det_old = det;
    M2 = M;
  }

  det = 1.;
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    alps::numeric::matrix<M_TYPE> G0(M[flavor].matrix());
    std::fill(G0.get_values().begin(), G0.get_values().end(), 0);

    const size_t Nv = M[flavor].matrix().num_rows();

    if (Nv==0) {
      M[flavor].matrix() = alps::numeric::matrix<M_TYPE>(0,0);
    } else {
      //construct G0 matrix
      for (size_t q=0; q<Nv; ++q) {
        for (size_t p=0; p<Nv; ++p) {
          G0(p, q) = mycast<M_TYPE>(green0_spline_for_M(flavor, p, q));
        }
      }
      for (size_t p=0; p<Nv; ++p) {
        G0(p, p) -= M[flavor].alpha_at(p);
      }
      if (verbose) {
        matrix_t MG = mygemm(M[flavor].matrix(),G0);
        double max_diff = is_identity(MG);
        if (max_diff>1E-8)
          std::cout<<" node= " << node << " step= " << step << " flavor " << flavor << " WARNING: roundoff errors in M*G, max_diff = " << max_diff << std::endl;
      }
      M[flavor].matrix() = alps::numeric::inverse(G0);
      det *= alps::numeric::determinant(G0);//the determinant is basically computed in alps::numeric::inverse(G0)
    }

    //reserve bit larger memory to avoid expensive memory reallocation
    const size_t size_reserved = num_cols(M[flavor].matrix())>0 ? std::min((size_t)max_order, (size_t)num_cols(1.5*M[flavor].matrix())) : 10;
    M[flavor].matrix().reserve(size_reserved, size_reserved);
  }

  //sign = boost::math::sign(det);
  //{
    //const std::vector<vertex_definition<M_TYPE> >& vd = Uijkl.get_vertices();
    //for (int iv=0; iv<itime_vertices.size(); ++iv)
      //sign *= boost::math::sign(-vd[itime_vertices[iv].type()].Uval());
  //}
  //std::cout << "sign, sign_det " << sign << " " << boost::math::sign(det) << std::endl;
  //assert(sign==boost::math::sign(det));

  //std::cout<<"debug : determinant computed by fast update = " << det_old << " determinant computed by matrix inversion = " << det << std::endl;

  if (verbose) {
    if (std::abs(det-det_old)/std::abs(det)>1E-8)
      std::cout<<" node= " << node << " step= " << step << " WARNING: roundoff errors : determinant computed by fast update = " << det_old << " determinant computed by matrix inversion = " << det << std::endl;
    for(unsigned int z=0;z<M2.size();++z){
      double max_diff=0;
      for(unsigned int j=0;j<num_cols(M2[z].matrix());++j){
        for(unsigned int i=0;i<num_rows(M2[z].matrix());++i){
          M_TYPE diff=M[z].matrix()(i,j)-M2[z].matrix()(i,j);
          if(std::abs(diff)>max_diff) max_diff=std::abs(diff);
        }
      }
      if(max_diff > 1.e-8)
        std::cout<<" node = " << node << " WARNING: roundoff errors in flavor: "<<z<<" max diff "<<max_diff<<std::endl;
    }
  }
}

