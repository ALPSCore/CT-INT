#include "../submatrix.hpp"


/***
 * Vertex insertion and removal
 */
template<typename T>
template<typename NVertexProb, typename R>
void SubmatrixUpdate<T>::vertex_insertion_removal_update(NVertexProb nv_prob, R& random) {
  pos_vertices_ins.resize(k_ins_max_);
  num_vertices_ins.resize(k_ins_max_);

  //add non-interacting vertices
  const int Nv0 = itime_vertices_.size();
  int vertex_begin = Nv0;
  for (int i_ins=0; i_ins<k_ins_max_; ++i_ins) {
    int Nv = nv_prob();
    assert(Nv>0);
    std::vector<itime_vertex> new_vertices;
    if (Nv==1) {
      new_vertices = generate_itime_vertices(*p_Uijkl_,random,beta_,Nv,all_type());
    } else {
      new_vertices = generate_itime_vertices(*p_Uijkl_,random,beta_,Nv,all_type());
      //throw std::runtime_error("Nv>2 is not implemented!");
    }

    for (int iv=0; iv<new_vertices.size(); ++iv) {
      new_vertices[iv].set_non_interacting();//this does not modify the spin state but just hide it.
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
  //const int try_ins = random()<0.5 ? 0 : 1;
  const int try_ins = 0; //DEBUG
  int i_ins = 0;
  for (int i_update=0; i_update<2*k_ins_max_; ++i_update) {
    if (i_update%2==try_ins) {
      std::cout << "trying: insertion " << i_ins <<  " nv " << num_vertices_ins[i_ins] << std::endl;
      assert(i_ins<pos_vertices_ins.size());
      insertion_step(random, pos_vertices_ins[i_ins], num_vertices_ins[i_ins]);
      ++i_ins;
    } else {
      std::cout << "trying: removal " << std::endl;
      removal_step(random, nv_prob()); //FOR DEBUG
    }
    for (int flavor=0; flavor<n_flavors(); ++flavor) {
      gamma_matrices_[flavor].sanity_check(invA_[flavor], spline_G0_);
    }
  }

  //update A^{-1}
  finalize_update();
};

template<typename T>
template<typename R>
void SubmatrixUpdate<T>::insertion_step(R& random, int vertex_begin, int num_vertices_ins) {
  //assert(vertex_begin<itime_vertices_.size());
  assert(vertex_begin+num_vertices_ins<=itime_vertices_.size());

  new_spins_work.resize(num_vertices_ins);
  pos_vertices_work.resize(num_vertices_ins);
  for (int iv=0; iv<num_vertices_ins; ++iv) {
    assert(iv+vertex_begin<itime_vertices_.size());
    new_spins_work[iv] = itime_vertices_[iv+vertex_begin].af_state();
    pos_vertices_work[iv] = iv+vertex_begin;
    //std::cout << "debug adding at time " << itime_vertices_[iv+vertex_begin].time() << std::endl;
    //std::cout << "debug adding with spin " << itime_vertices_[iv+vertex_begin].af_state() << std::endl;
  }

  const T det_rat_A = try_spin_flip(pos_vertices_work, new_spins_work);
  //const T det_rat_M = coeff_det*det_rat_A;

  //actually flip spins to compute proposal rates
  //for (int iv=0; iv<num_vertices_ins; ++iv) {
    //itime_vertices_[iv + vertex_begin].set_interacting();
  //}

  /****** TO DO compute prod_U *******/

  //if (std::abs(det_rat_A)>random()) {
  if (random()>0.5) {//DEBUG
    std::cout << "accepted " << std::endl;
    det_A_ *= det_rat_A;
    perform_spin_flip(pos_vertices_work, new_spins_work);
  } else {
    std::cout << "rejected " << std::endl;
    reject_spin_flip();
    //for (int iv=0; iv<num_vertices_ins; ++iv) {
      //itime_vertices_[iv + vertex_begin].set_non_interacting();
    //}
  }
};

template<typename T>
template<typename R>
void SubmatrixUpdate<T>::removal_step(R& random, int nv_rem) {
  const int Nv = itime_vertices_.size();
  std::vector<int> pos_int_vertices(Nv);
  int num_int = 0;
  for (int iv=0; iv<Nv; ++iv) {
    if (!itime_vertices_[iv].is_non_interacting()) {
      pos_int_vertices[num_int] = iv;
      ++num_int;
    }
  }

  //std::cout << "debug nv_rem " << nv_rem << std::endl;

  if (num_int<nv_rem) return;//there are not an enough number of vertices to removed.

  std::vector<int> rand_pos = pickup_a_few_numbers(num_int, nv_rem, random);

  //for (int iv=0; iv<num_int; ++iv) {
    //std::cout << "int " << iv << " " << pos_int_vertices[iv] << std::endl;
  //}

  std::vector<int> pos_vertices_remove(nv_rem);
  for (int iv=0; iv<nv_rem; ++iv) {
    pos_vertices_remove[iv] = pos_int_vertices[rand_pos[iv]];
    //std::cout << "debug removing at time " << itime_vertices_[pos_vertices_remove[iv]].time() << std::endl;
    assert(pos_vertices_remove[iv]<Nv);
    assert(!itime_vertices_[pos_vertices_remove[iv]].is_non_interacting());
  }
  std::vector<int> new_spins_remove(nv_rem, NON_INT_SPIN_STATE);

  const T det_rat_A = try_spin_flip(pos_vertices_remove, new_spins_remove);

  //if (std::abs(det_rat_A)>random()) {
  if (random()>0.5) {//DEBUG
    std::cout << "accepted " << std::endl;
    det_A_ *= det_rat_A;
    perform_spin_flip(pos_vertices_remove, new_spins_remove);
  } else {
    std::cout << "rejected " << std::endl;
    reject_spin_flip();
  }
};

