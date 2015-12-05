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
void InteractionExpansion::removal_insertion_update(void) {
  const int nv_updated = update_prop.gen_Nv(random.engine());

  if (nv_updated==1) {
    removal_insertion_single_vertex_update();
  } else if (nv_updated==2) {
    removal_insertion_double_vertex_update();
  } else {
    throw std::runtime_error("Not implemented");
  }

  for(spin_t flavor=0; flavor<n_flavors; ++flavor) {
    vertex_histograms[flavor]->count(num_rows(M[flavor].matrix()));
  }
}

/*
void InteractionExpansion::removal_insertion_update(void)
{
  const int nv_updated = update_prop.gen_Nv(random.engine());

  const int pert_order= itime_vertices.size();   //current order of perturbation series
  double metropolis_weight=0.;
  double det_rat=0;
  if(random()<0.5){  //trying to ADD vertex
    M.sanity_check(itime_vertices);
    if(pert_order+nv_updated>max_order)
      return; //we have already reached the highest perturbation order
    if (nv_updated>=2) {
      add_helper.op = non_density_type_in_window(beta*random(),
                                                 std::min(beta,window_width+window_dist(random.engine())),
                                                 beta);
    }
    std::vector<itime_vertex> new_vertices;
    std::pair<int,int> v_pair;
    if (nv_updated==1) {
      if (single_vertex_update_non_density_type) {
        new_vertices = generate_itime_vertices(Uijkl,random,beta,nv_updated,all_type());
      } else {
        new_vertices = generate_itime_vertices(Uijkl,random,beta,nv_updated,density_type());
      }
    } else if (nv_updated==2) {
      v_pair = mv_update_valid_pair[mv_update_valid_pair.size()*random()];
      new_vertices = generate_valid_vertex_pair2(Uijkl,v_pair,random,beta,symm_exp_dist);
    } else {
      throw std::runtime_error("Not implemented!");
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
    boost::tie(metropolis_weight,det_rat)=try_add(add_helper, nv_updated, new_vertices);

    double p_ins, p_rem;
    if (nv_updated==1) {
      if (single_vertex_update_non_density_type) {
        p_ins = 1.0 / (beta * Uijkl.n_vertex_type());
        p_rem = 1.0 / itime_vertices_new.size();
      } else {
        p_ins = 1.0 / (beta * Uijkl.num_density_vertex_type());
        p_rem = 1.0 / std::count_if(itime_vertices_new.begin(),itime_vertices_new.end(),density_type());
      }
    } else if (nv_updated==2) {
      p_ins = 1.;
      const double dtau = mymod(new_vertices[1].time()-new_vertices[0].time(), beta);
      int n_vpair;
      double F;
      pick_up_valid_vertex_pair2(itime_vertices_new,
                                 std::make_pair(new_vertices[0].type(),new_vertices[1].type()),
                                 beta, symm_exp_dist, random, n_vpair, F);
      p_rem = (beta*beta)*symm_exp_dist.coeff_X(dtau)/F;
      if (n_vpair==0)
        throw std::logic_error("v_pair must be larger than 0.");
    } else {
      throw std::runtime_error("Not implemented!");
    }
    metropolis_weight *= p_rem/p_ins;

    if (nv_updated>=2) {
      statistics_ins.add_sample(compute_spread(new_vertices,beta), std::min(fabs(metropolis_weight),1.0), nv_updated-2);
    }
    if (nv_updated==2) {
      statistics_dv_ins.add_sample(mymod(new_vertices[1].time()-new_vertices[0].time(), beta),
                                   std::min(fabs(metropolis_weight), 1.0),
                                   std::distance(mv_update_valid_pair.begin(), std::find(mv_update_valid_pair.begin(), mv_update_valid_pair.end(), v_pair))
      );
    }
    if(fabs(metropolis_weight)> random()){
      perform_add(add_helper,nv_updated);
      sign*=boost::math::sign(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      simple_statistics_ins.accepted(nv_updated-1);
      if(fabs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);
    }else{
      reject_add(add_helper,nv_updated);
      M.sanity_check(itime_vertices);
      simple_statistics_ins.rejected(nv_updated-1);
    }
#ifndef NDEBUG
    sanity_check();
#endif
  }else{ // try to REMOVE a vertex
    //std::cout << "try rem " << std::endl;
    M.sanity_check(itime_vertices);
    if(pert_order < nv_updated) {
      return;
    }
    if (nv_updated>=2) {
      remove_helper.op = non_density_type_in_window(beta*random(),
              std::min(beta,window_width+window_dist(random.engine())),
              beta);
    }

    //choose vertices to be removed
    std::vector<int> vertices_nr;
    int n_vpair;
    double F;
    double p_ins, p_rem;
    if (nv_updated==1) {
      if (single_vertex_update_non_density_type) {
        vertices_nr = pick_up_itime_vertices(itime_vertices, random, nv_updated, all_type());
      } else {
        vertices_nr = pick_up_itime_vertices(itime_vertices, random, nv_updated, density_type());
      }
    } else if (nv_updated==2) {
      std::pair<int,int> v_pair = mv_update_valid_pair[mv_update_valid_pair.size()*random()];
      std::pair<int,int> r = pick_up_valid_vertex_pair2(itime_vertices, v_pair, beta, symm_exp_dist, random, n_vpair, F);
      //std::cout << "debug " << itime_vertices.size() << " " << r.first << "  " << r.second << std::endl;
      if (n_vpair>0) {
        vertices_nr.resize(2);
        vertices_nr[0]=r.first;
        vertices_nr[1]=r.second;
      } else {
        vertices_nr.resize(0);
      }
    } else {
      throw std::runtime_error("Not implemented!");
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

    if (nv_updated==1) {
      if (single_vertex_update_non_density_type) {
        p_ins = 1.0/(beta*Uijkl.n_vertex_type());
        p_rem = 1.0/pert_order;
      } else {
        p_ins = 1.0/(beta*Uijkl.num_density_vertex_type());
        p_rem = 1.0/std::count_if(itime_vertices.begin(), itime_vertices.end(), density_type());
      }
    } else if (nv_updated==2) {
      assert(vertices_to_be_removed[1].type()>vertices_to_be_removed[0].type());
      p_ins = 1.0;
      p_rem = (beta*beta)*
        symm_exp_dist.coeff_X(
          mymod(vertices_to_be_removed[1].time()-vertices_to_be_removed[0].time(), beta)
        )/F;
    }

    update_prop.generated_valid_update(nv_updated);

    boost::tie(metropolis_weight,det_rat)=try_remove(vertices_nr, remove_helper); //get the determinant ratio. don't perform fastupdate yet
    metropolis_weight *= p_ins/p_rem;
    if (nv_updated>=2) {
      statistics_rem.add_sample(compute_spread(vertices_to_be_removed, beta), std::min(fabs(metropolis_weight), 1.0),
                                nv_updated - 2);
    }
    if (nv_updated==2) {
      int idx = std::distance(mv_update_valid_pair.begin(),
                              std::find(
                                  mv_update_valid_pair.begin(),
                                  mv_update_valid_pair.end(),
                                  std::make_pair(vertices_to_be_removed[0].type(),vertices_to_be_removed[1].type())
                              )
      );
      statistics_dv_rem.add_sample(mymod(vertices_to_be_removed[1].time()-vertices_to_be_removed[0].time(), beta), std::min(fabs(metropolis_weight), 1.0), idx);
    }
    if(fabs(metropolis_weight)> random()){ //do the actual update
      perform_remove(vertices_nr, remove_helper);
      sign*=boost::math::sign(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      simple_statistics_rem.accepted(nv_updated-1);

      if(fabs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);

      //std::cout << " qn = " << is_quantum_number_conserved(itime_vertices) << std::endl;
      //if(fabs(metropolis_weight)<1E-3 || (step>= 2558981 && step<=2559000) )  {
        //std::cout << "doing sanity check " << metropolis_weight << " node " << node << std::endl;
        //reset_perturbation_series(true);
        //std::cout << "sanity check done " << metropolis_weight << " node " << node << std::endl;
      //}

    }else{
      //measurements["VertexRemoval"]<<0.;
      reject_remove(remove_helper);
      M.sanity_check(itime_vertices);
      simple_statistics_rem.rejected(nv_updated-1);
    }
#ifndef NDEBUG
    sanity_check();
#endif
  }//end REMOVE
  weight=metropolis_weight;
  for(spin_t flavor=0; flavor<n_flavors; ++flavor) {
    vertex_histograms[flavor]->count(num_rows(M[flavor].matrix()));
  }
}
*/

void InteractionExpansion::removal_insertion_single_vertex_update(void)
{
  const int nv_updated = 1;

  const int pert_order= itime_vertices.size();   //current order of perturbation series
  double metropolis_weight=0.;
  double det_rat=0;
  if(random()<0.5){  //trying to ADD vertex
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
    boost::tie(metropolis_weight,det_rat)=try_add(add_helper, nv_updated, new_vertices);

    double p_ins, p_rem;
    if (single_vertex_update_non_density_type) {
      p_ins = 1.0 / (beta * Uijkl.n_vertex_type());
      p_rem = 1.0 / itime_vertices_new.size();
    } else {
      p_ins = 1.0 / (beta * Uijkl.num_density_vertex_type());
      p_rem = 1.0 / std::count_if(itime_vertices_new.begin(),itime_vertices_new.end(),density_type());
    }
    metropolis_weight *= p_rem/p_ins;

    if(fabs(metropolis_weight)> random()){
      perform_add(add_helper,nv_updated);
      sign*=boost::math::sign(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      simple_statistics_ins.accepted(nv_updated-1);
      if(fabs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);
    }else{
      reject_add(add_helper,nv_updated);
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

    boost::tie(metropolis_weight,det_rat)=try_remove(vertices_nr, remove_helper); //get the determinant ratio. don't perform fastupdate yet
    metropolis_weight *= p_ins/p_rem;
    if(fabs(metropolis_weight)> random()){ //do the actual update
      perform_remove(vertices_nr, remove_helper);
      sign*=boost::math::sign(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      simple_statistics_rem.accepted(nv_updated-1);
      if(fabs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);
    }else{
      reject_remove(remove_helper);
      M.sanity_check(itime_vertices);
      simple_statistics_rem.rejected(nv_updated-1);
    }
#ifndef NDEBUG
    sanity_check();
#endif
  }//end REMOVE
  weight=metropolis_weight;
}

void InteractionExpansion::removal_insertion_double_vertex_update(void)
{
  const int nv_updated = 2;

  if (mv_update_valid_pair.size()==0)
    return;

  const int pert_order= itime_vertices.size();   //current order of perturbation series
  double metropolis_weight=0.;
  double det_rat=0;
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
    boost::tie(metropolis_weight,det_rat)=try_add(add_helper, nv_updated, new_vertices);

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

    statistics_ins.add_sample(compute_spread(new_vertices,beta), std::min(fabs(metropolis_weight),1.0), nv_updated-2);
    statistics_dv_ins.add_sample(mymod(new_vertices[1].time()-new_vertices[0].time(), beta),
                                 std::min(fabs(metropolis_weight), 1.0),
                                 std::distance(mv_update_valid_pair.begin(), std::find(mv_update_valid_pair.begin(), mv_update_valid_pair.end(), v_pair))
    );
    if(fabs(metropolis_weight)> random()){
      perform_add(add_helper,nv_updated);
      sign*=boost::math::sign(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      simple_statistics_ins.accepted(nv_updated-1);
      if(fabs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);
    }else{
      reject_add(add_helper,nv_updated);
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

    boost::tie(metropolis_weight,det_rat)=try_remove(vertices_nr, remove_helper); //get the determinant ratio. don't perform fastupdate yet
    metropolis_weight /= (beta*beta)*
                         symm_exp_dist.coeff_X(
                             mymod(vertices_to_be_removed[1].time()-vertices_to_be_removed[0].time(), beta)
                         )/F;
    statistics_rem.add_sample(compute_spread(vertices_to_be_removed, beta), std::min(fabs(metropolis_weight), 1.0),
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
        std::min(fabs(metropolis_weight), 1.0),
        idx);
    if(fabs(metropolis_weight)> random()){ //do the actual update
      perform_remove(vertices_nr, remove_helper);
      sign*=boost::math::sign(metropolis_weight);
      det*=det_rat;
      M.sanity_check(itime_vertices);
      simple_statistics_rem.accepted(nv_updated-1);
      if(fabs(metropolis_weight)<1E-5)
        reset_perturbation_series(false);
    }else{
      reject_remove(remove_helper);
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
void InteractionExpansion::shift_update(void) {
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

  double metropolis_weight = try_shift(iv, new_time);

  statistics_shift.add_sample(std::min(diff_time,beta-diff_time), std::min(fabs(metropolis_weight),1.0), 0);
  if(fabs(metropolis_weight)> random()){ //do the actual update
    perform_shift(iv);
    sign*=boost::math::sign(metropolis_weight);
    det*=metropolis_weight;
    if(fabs(metropolis_weight)<1E-5)
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

//Update alpha for non-density-type vertices
void InteractionExpansion::alpha_update(void) {
  if (itime_vertices.size()==0 || alpha_scale_max==1.0)
    return;

  std::cout << "current alpha_scale " << alpha_scale << std::endl;
  double new_alpha_scale = std::exp(random()*(std::log(alpha_scale_max)-std::log(alpha_scale_min))+std::log(alpha_scale_min));
  const double ratio = 1.05;
  //const double new_alpha_scale = random()<0.5 ? std::min(alpha_scale*ratio, alpha_scale_max) : std::max(alpha_scale/ratio, alpha_scale_min);
  //std::cout <<  "random " << 2*(static_cast<int>(random())-1)*boost::random::exponential_distribution<>(1/std::log(ratio))(random.engine()) << std::endl;
  //double new_alpha_scale = std::exp(
      //2*(static_cast<int>(random())-1)*boost::random::exponential_distribution<>(1/std::log(ratio))(random.engine())+std::log(alpha_scale)
  //);
  //double new_alpha_scale = random()<0.5 ? std::min(alpha_scale*ratio, alpha_scale_max) : std::max(alpha_scale/ratio, alpha_scale_min);
  if (new_alpha_scale>alpha_scale_max || new_alpha_scale <alpha_scale_min)
    return;

  const double det_old = det;
  const double sign_old = sign;
  const big_inverse_m_matrix M_old(M);

  M.set_alpha_scale(new_alpha_scale);
  boost::timer::cpu_timer timer;
  reset_perturbation_series(false);
  const double metropolis_weight = (det/det_old) * std::exp(-flat_histogram_alpha.log_dos_estimated(std::log(new_alpha_scale))+flat_histogram_alpha.log_dos_estimated(std::log(alpha_scale)));
  std::cout << "debug " << alpha_scale << " " << new_alpha_scale
    << " " << flat_histogram_alpha.log_dos_estimated(std::log(alpha_scale))
    << " " << flat_histogram_alpha.log_dos_estimated(std::log(new_alpha_scale))
    << " " << flat_histogram_alpha.to_idx(std::log(alpha_scale))
    << " " << flat_histogram_alpha.to_idx(std::log(new_alpha_scale)) << std::endl;
  //std::cout << " timing_alpha_update " << timer.elapsed().wall*1E-6 << " msec. pert_order = " << itime_vertices.size() << std::endl;

  std::cout << "debug " << metropolis_weight << " new_alpha_scale " << new_alpha_scale << " qn " << is_quantum_number_conserved(itime_vertices) << std::endl;
  if(fabs(metropolis_weight)> random()){ //do the actual update
    alpha_scale = new_alpha_scale;
    //std::cout << "node " << node << " step " << step << " alpha_update accepted new alpha_scale = " << alpha_scale << " " << fabs(metropolis_weight) << std::endl;
  }else {
    det = det_old;
    M = M_old;
    sign = sign_old;
    //std::cout << "node " << node << " step " << step << " alpha_update rejected new alpha_scale = " << alpha_scale << std::endl;
    //std::cout << "step " << step << " alpha_update rejected " << std::endl;
  }
}

///We recreate M from scratch to avoid roundoff error.
//This is done by constructing a matrix consisting of bare Green's function and invert it.
void InteractionExpansion::reset_perturbation_series(bool verbose)
{
  typedef alps::numeric::matrix<GTYPE> matrix_t;
  big_inverse_m_matrix M2;
  double det_old;
  if (verbose) {
    det_old = det;
    M2 = M;
  }

  det = 1.;
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
  }

  sign = boost::math::sign(det);
  {
    const std::vector<vertex_definition<GTYPE> >& vd = Uijkl.get_vertices();
    for (int iv=0; iv<itime_vertices.size(); ++iv)
      sign *= boost::math::sign(-vd[itime_vertices[iv].type()].Uval());
  }
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
          double diff=M[z].matrix()(i,j)-M2[z].matrix()(i,j);
          if(std::abs(diff)>max_diff) max_diff=std::abs(diff);
        }
      }
      if(max_diff > 1.e-8)
        std::cout<<" node = " << node << " WARNING: roundoff errors in flavor: "<<z<<" max diff "<<max_diff<<std::endl;
    }
  }
}

