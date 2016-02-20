#include "../submatrix.hpp"


/***
 * Vertex insertion and removal
 */
/*
template<typename T>
template<typename MANAGER, typename R>
T SubmatrixUpdate<T>::vertex_insertion_removal_update(MANAGER& manager, R& random) {

  int num_ins_try = 0;
  std::vector<bool> try_ins(2*k_ins_max_, false);
  for (int i_update=0; i_update<2*k_ins_max_; ++i_update) {
    if (random()<0.5) {
      try_ins[i_update] = true;
      ++num_ins_try;
    }
  }

  pos_vertices_ins.resize(num_ins_try);
  num_vertices_ins.resize(num_ins_try);

  T weight_rat = 1.0;

  std::vector<typename MANAGER::insertion_removal_info_type> ins_info(num_ins_try);

  //add non-interacting vertices
  const int Nv0 = itime_vertices_.size();
  int vertex_begin = Nv0;
  for (int i_ins=0; i_ins<num_ins_try; ++i_ins) {

    std::vector<itime_vertex> new_vertices;
    boost::tie(new_vertices,ins_info[i_ins]) = manager.gen_itime_vertices_insertion(*p_Uijkl_, random);
    const int Nv = new_vertices.size();

    for (int iv=0; iv<new_vertices.size(); ++iv) {
      new_vertices[iv].set_non_interacting();//this does not modify the spin state but just hide it.
      new_vertices[iv].set_unique_id(gen_new_vertex_id());
      itime_vertices_.push_back(new_vertices[iv]);
    }
    pos_vertices_ins[i_ins] = vertex_begin;
    num_vertices_ins[i_ins] = new_vertices.size();
    vertex_begin += new_vertices.size();
  }
  assert(pos_vertices_ins.size()==num_vertices_ins.size());

  //set starting point
  init_update(Nv0);

  //perform actual updates
  int i_ins = 0;
  for (int i_update=0; i_update<2*k_ins_max_; ++i_update) {
    T rtmp;
    if (try_ins[i_update]) {
      assert(i_ins<pos_vertices_ins.size());
      rtmp = insertion_step(manager, random, pos_vertices_ins[i_ins], num_vertices_ins[i_ins]);
      manager.feedback_insertion(ins_info[i_ins], rtmp!=1.0);
      ++i_ins;
    } else {
      rtmp = removal_step(manager, random);
    }
    weight_rat *= rtmp;
#ifndef NDEBUG
    for (int flavor=0; flavor<n_flavors(); ++flavor) {
      gamma_matrices_[flavor].sanity_check(invA_[flavor], spline_G0_);
    }
#endif
  }
  assert(i_ins==num_ins_try);

  //update A^{-1}
  finalize_update();

  return weight_rat;
};

template<typename T>
template<typename MANAGER, typename R>
T SubmatrixUpdate<T>::insertion_step(MANAGER& manager, R& random, int vertex_begin, int num_vertices_ins) {
  assert(vertex_begin+num_vertices_ins<=itime_vertices_.size());

  if (num_vertices_ins==0) {
    return 1.0;
  }

  T det_rat_A, f_rat, U_rat;

  new_spins_work.resize(num_vertices_ins);
  pos_vertices_work.resize(num_vertices_ins);
  for (int iv=0; iv<num_vertices_ins; ++iv) {
    assert(iv+vertex_begin<itime_vertices_.size());
    assert(itime_vertices_[iv+vertex_begin].is_non_interacting());
    new_spins_work[iv] = itime_vertices_[iv+vertex_begin].af_state();
    pos_vertices_work[iv] = iv+vertex_begin;
  }

  boost::tie(det_rat_A,f_rat,U_rat) = try_spin_flip(pos_vertices_work, new_spins_work);
  const double acc_corr = manager.acc_rate_corr_insertion(itime_vertices_, pos_vertices_work, random);

  T prob = det_rat_A*f_rat*U_rat*acc_corr;

  if (std::abs(prob)>random()) {
    //std::cout << "accepted " << std::endl;
    perform_spin_flip(pos_vertices_work, new_spins_work);
    return det_rat_A*f_rat*U_rat;
  } else {
    //std::cout << "rejected " << std::endl;
    reject_spin_flip();
    return 1.0;
  }
}

template<typename T>
template<typename MANAGER, typename R>
T SubmatrixUpdate<T>::removal_step(MANAGER& manager, R& random) {
  const int Nv = itime_vertices_.size();
  std::vector<int> pos_vertices_remove;
  const double acc_corr = manager.pick_up_vertices_to_be_removed(itime_vertices_, random, pos_vertices_remove);
  const int nv_rem = pos_vertices_remove.size();

  //no vertices to be removed
  if (nv_rem==0) {
    return 1.0;
  }

#ifndef NDEBUG
  for (int iv=0; iv<pos_vertices_remove.size(); ++iv) {
    assert(!itime_vertices_[pos_vertices_remove[iv]].is_non_interacting());
  }
#endif

  std::vector<int> new_spins_remove(nv_rem, NON_INT_SPIN_STATE);

  T det_rat_A, f_rat, U_rat;
  boost::tie(det_rat_A,f_rat,U_rat) = try_spin_flip(pos_vertices_remove, new_spins_remove);

  T prob = det_rat_A*f_rat*U_rat*acc_corr;

  if (std::abs(prob)>random()) {
    //std::cout << "accepted " << std::endl;
    perform_spin_flip(pos_vertices_remove, new_spins_remove);
    return det_rat_A*f_rat*U_rat;
  } else {
    //std::cout << "rejected " << std::endl;
    reject_spin_flip();
    return 1.0;
  }
}
 */
