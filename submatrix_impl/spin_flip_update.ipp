#include "../submatrix.hpp"


/***
 * Spin flip update
 */
template<typename T>
template<typename R>
T SubmatrixUpdate<T>::spin_flip_update(R& random) {
  //pos_vertices_ins.resize(k_ins_max_);
  //num_vertices_ins.resize(k_ins_max_);

  const int Nv0 = itime_vertices_.size();
  const int nv_flip = std::min(2*k_ins_max_, Nv0);

  T weight_rat = 1.0;

  std::vector<int> pos_vertices_flip = pickup_a_few_numbers(Nv0, nv_flip, random);
  assert(pos_vertices_flip.size()==nv_flip);

  //set starting point (do not add non-interacting vertices)
  init_update(Nv0);

  //perform actual updates
  for (int i_update=0; i_update<nv_flip; ++i_update) {
    weight_rat *= spin_flip_step(random, pos_vertices_flip[i_update]);
  }

  //update A^{-1}
  finalize_update();

  return weight_rat;
};

template<typename T>
template<typename R>
T SubmatrixUpdate<T>::spin_flip_step(R& random, int pos_vertex) {
  T det_rat_A, f_rat, U_rat;

  std::vector<int> pos_vertices_tmp(1), new_spins_tmp(1);
  pos_vertices_tmp[0] = pos_vertex;
  const int num_af_states = p_Uijkl_->get_vertex(itime_vertices_[pos_vertex].type()).num_af_states();

  if (num_af_states<=1) {
    return 1.0;
  }

  new_spins_tmp[0] = (int) (random()*num_af_states);
  while (new_spins_tmp[0] == itime_vertices_[pos_vertex].af_state()) {
    new_spins_tmp[0] = (int) (random()*num_af_states);
  }

  boost::tie(det_rat_A,f_rat,U_rat) = try_spin_flip(pos_vertices_tmp, new_spins_tmp);
  const T prob = det_rat_A*f_rat*U_rat;

  if (std::abs(prob)>random()) {
    //std::cout << "accepted " << std::endl;
    perform_spin_flip(pos_vertices_tmp, new_spins_tmp);
    return prob;
  } else {
    //std::cout << "rejected " << std::endl;
    reject_spin_flip();
    return 1.0;
  }
}
