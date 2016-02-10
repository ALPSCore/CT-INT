#include "../submatrix.hpp"

template<typename T, typename SPLINE_G0_TYPE>
T eval_Gij(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0, int row_A, int col_A) {
  static alps::numeric::matrix<T> G0, invA_G0(1,1);

  T alpha_col = invA.alpha_at(col_A);
  if (alpha_col!=ALPHA_NON_INT) {
    //use Eq. (A3)
    const T fj = eval_f(alpha_col);
    return row_A==col_A ? (fj*invA.matrix()(row_A,col_A)-1.0)/(fj-1.0)
                        : (fj*invA.matrix()(row_A,col_A))/(fj-1.0);
  } else {
    //use Eq. (A4)
    const int Nv = invA.matrix().num_rows();
    assert (invA.creators().size()==Nv);
    assert (invA.annihilators().size()==Nv);
    G0.resize_values_not_retained(Nv, 1);
    for (int iv=0; iv<Nv; ++iv) {
      G0(iv,0) = spline_G0(invA.annihilators()[iv], invA.creators()[col_A]);
    }
    alps::numeric::submatrix_view<T> invA_view(invA.matrix(), row_A, 0, 1, Nv);
    mygemm((T) 1.0, invA_view, G0, (T) 0.0, invA_G0);
    return invA_G0(0,0);
  }
}

template<typename T>
void SubmatrixUpdate<T>::sanity_check() const {
#ifndef NDEBUG
  //check gamma^{-1}
  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    gamma_matrices_[flavor].sanity_check(invA_[flavor], spline_G0_);
  }

  //check A^{-1}
  if (state==READY_FOR_UPDATE) {
    invA_.sanity_check(spline_G0_, p_Uijkl_, itime_vertices_);
  }
#endif
}

template<typename T>
void SubmatrixUpdate<T>::recompute(bool check_error) {
  if (state==READY_FOR_UPDATE) {
    invA_.recompute(spline_G0_, check_error);
  }
}

template<typename T>
SubmatrixUpdate<T>::SubmatrixUpdate(int k_ins_max, int n_flavors, SPLINE_G0_TYPE spline_G0, general_U_matrix<T>* p_Uijkl, double beta) :
    k_ins_max_(k_ins_max),
    spline_G0_(spline_G0),
    p_Uijkl_(p_Uijkl),
    beta_(beta),
    coeff_det(n_flavors%2==0 ? 1 : -1),
    alpha_scale_(1),
    invA_(n_flavors),
    det_A_(1.0),
    state(READY_FOR_UPDATE),
    gamma_matrices_(n_flavors),
    ops_rem(n_flavors),
    ops_ins(n_flavors),
    ops_replace(n_flavors)
{
}

template<typename T>
SubmatrixUpdate<T>::SubmatrixUpdate(int k_ins_max, int n_flavors, SPLINE_G0_TYPE spline_G0, general_U_matrix<T>* p_Uijkl, double beta,
                                    const itime_vertex_container& itime_vertices_init) :
    k_ins_max_(k_ins_max),
    spline_G0_(spline_G0),
    p_Uijkl_(p_Uijkl),
    beta_(beta),
    coeff_det(n_flavors%2==0 ? 1 : -1),
    alpha_scale_(1),
    invA_(n_flavors),
    det_A_(1.0),
    state(READY_FOR_UPDATE),
    gamma_matrices_(n_flavors),
    ops_rem(n_flavors),
    ops_ins(n_flavors),
    ops_replace(n_flavors)
{
  if (itime_vertices_init.size()!=0) {
    invA_.add_interacting_vertices(p_Uijkl, spline_G0, itime_vertices_init, 0);
    itime_vertices_ = itime_vertices_init;
  }
}

template<typename T>
void SubmatrixUpdate<T>::init_update(int begin_index) {
  assert(state==READY_FOR_UPDATE);

  //set starting point
  itime_vertices0_ = itime_vertices_;// tilde C^0 in Eq. (20)

  for (int iv=0; iv<begin_index; ++iv) {
    assert(!itime_vertices0_[iv].is_non_interacting());
  }
  for (int iv=begin_index; iv<itime_vertices0_.size(); ++iv) {
    assert(itime_vertices0_[iv].is_non_interacting());
  }

  //extend the size of A^{-1} if we have non-interacting vertices
  if (begin_index<itime_vertices_.size()) {
    invA_.add_non_interacting_vertices(p_Uijkl_, spline_G0_, itime_vertices0_, begin_index);
  }

  state = TRYING_SPIN_FLIP;
}

template<typename T>
void SubmatrixUpdate<T>::finalize_update() {
  assert(state==TRYING_SPIN_FLIP);

  invA_.update_matrix(gamma_matrices_, spline_G0_);

  //remove cols and rows corresponding to non-interacting vertices
  std::vector<double> time_removed;
  itime_vertex_container new_itime_vertices;
  for (int iv=0; iv<itime_vertices_.size(); ++iv) {
    if (itime_vertices_[iv].is_non_interacting()) {
      time_removed.push_back(itime_vertices_[iv].time());
    } else {
      new_itime_vertices.push_back(itime_vertices_[iv]);
    }
  }
  invA_.remove_rows_cols(time_removed);
  std::swap(itime_vertices_, new_itime_vertices);

  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    gamma_matrices_[flavor].clear();
  }

  state = READY_FOR_UPDATE;
}

//returns the determinant ratio of A.
template<typename T>
T SubmatrixUpdate<T>::try_spin_flip(const std::vector<int>& pos, const std::vector<int>& new_spins) {
  assert(state==TRYING_SPIN_FLIP);

  const int n_pos = pos.size();

  ops_rem.resize(n_flavors());
  ops_ins.resize(n_flavors());
  ops_replace.resize(n_flavors());
  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    ops_rem[flavor].resize(0);
    ops_ins[flavor].resize(0);
    ops_replace[flavor].resize(0);
  }

#ifndef NDEBUG
  assert(pos.size()==new_spins.size());
  for (int i=0; i<n_pos; ++i) {
    assert (pos[i]<itime_vertices_.size());
  }
  std::vector<bool> visited(itime_vertices_.size(), false);
  for (int i=0; i<n_pos; ++i) {
    //std::cout << "i_pos " << i << " " << pos[i] << std::endl;
    assert(visited[pos[i]]==false);
    visited[pos[i]] = true;
  }
#endif

  //figure out which rows/cols are to be updated
  for (int ipos=0; ipos<pos.size(); ++ipos) {
    const int iv = pos[ipos];
    const int new_spin = new_spins[ipos];
    assert(iv<itime_vertices0_.size());
    const itime_vertex& v = itime_vertices0_[iv];
    const vertex_definition<T>& vdef = p_Uijkl_->get_vertex(v.type());
    for (int rank=0; rank<vdef.rank(); ++rank) {
      const int flavor_rank = vdef.flavors()[rank];
      const operator_time op_t = operator_time(v.time(),-rank);
      const T alpha0 = itime_vertices0_[iv].is_non_interacting()
                       ? ALPHA_NON_INT
                       : vdef.get_alpha(itime_vertices0_[iv].af_state(),rank);
      const T alpha_current = itime_vertices_[iv].is_non_interacting()
                       ? ALPHA_NON_INT
                       : vdef.get_alpha(itime_vertices_[iv].af_state(),rank);
      const T alpha_new = new_spin == NON_INT_SPIN_STATE
                       ? ALPHA_NON_INT
                       : vdef.get_alpha(new_spin, rank);
      const int pos_in_A = invA_[flavor_rank].find_row_col(op_t.time());
      assert(pos_in_A>=0);
      //std::cout << "ipos, rank " << ipos << " " << rank << " flavor " << flavor_rank << " " << pos_in_A <<std::endl;
      OperatorToBeUpdated<T> elem(op_t, pos_in_A, alpha0, alpha_current, alpha_new);

      if (new_spin==NON_INT_SPIN_STATE) {
        //if(itime_vertices0_[iv].is_non_interacting()) {
        if(alpha0==ALPHA_NON_INT) {
          ops_rem[flavor_rank].push_back(elem);
          //std::cout << "debug A " << std::endl;
        } else {
          ops_ins[flavor_rank].push_back(elem);
          //std::cout << "debug B " << std::endl;
        }
      } else {
        //new spin state is physical (interacting)
        if(alpha0==ALPHA_NON_INT) {
          ops_ins[flavor_rank].push_back(elem);
          //std::cout << "debug C " << std::endl;
        } else {
          if (alpha0 != alpha_new) {
            ops_replace[flavor_rank].push_back(elem);
            //std::cout << "debug D " << std::endl;
          }
        }
      }
    }
  }

  //compute determinant ratio of A (not M):
  T det_rat_A = (T) 1.0;
  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    //std::cout << "debug ";
    //std::cout << ops_rem.size();
    //std::cout << ops_ins.size();
    //std::cout << ops_replace.size();
    //std::cout << std::endl;

    if (ops_rem.size()==0 && ops_ins.size()==0 && ops_replace.size()==0) {
      continue;
    } else if (ops_rem[flavor].size()==0 && ops_ins[flavor].size()>0 && ops_replace[flavor].size()==0) {
        det_rat_A *= gamma_matrices_[flavor].try_add(invA_[flavor], spline_G0_, ops_ins[flavor]);
    } else if (ops_rem[flavor].size()>0 && ops_ins[flavor].size()==0 && ops_replace[flavor].size()==0) {
      det_rat_A *= gamma_matrices_[flavor].try_remove(invA_[flavor], spline_G0_, ops_rem[flavor]);
    } else if (ops_rem[flavor].size()==0 && ops_ins[flavor].size()==0 && ops_replace[flavor].size()>0) {
      throw std::runtime_error("Not implemented: try_replace");
      //det_rat_A *= gamma_matrices_[flavor].try_replace(p_Uijkl_, invA_[flavor], ops_ins[flavor]);
    } else {
      throw std::runtime_error("Not implemented: try_***");
    }
  }
  return det_rat_A;
}

template<typename T>
void SubmatrixUpdate<T>::perform_spin_flip(const std::vector<int>& pos, const std::vector<int>& new_spins) {
  assert(state==TRYING_SPIN_FLIP);

  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    if (ops_rem.size()==0 && ops_ins.size()==0 && ops_replace.size()==0) {
      continue;
    } else if (ops_rem[flavor].size()==0 && ops_ins[flavor].size()>0 && ops_replace[flavor].size()==0) {
      gamma_matrices_[flavor].perform_add();
    } else if (ops_rem[flavor].size()>0 && ops_ins[flavor].size()==0 && ops_replace[flavor].size()==0) {
      gamma_matrices_[flavor].perform_remove();
    } else if (ops_rem[flavor].size()==0 && ops_ins[flavor].size()==0 && ops_replace[flavor].size()>0) {
      throw std::runtime_error("Not implemented");
      //gamma_matrices_[flavor].perform_replace();
    } else {
      throw std::runtime_error("Not implemented");
    }
  }

  for (int i_pos=0; i_pos<pos.size(); ++i_pos) {
    if (new_spins[i_pos]==NON_INT_SPIN_STATE) {
      itime_vertices_[pos[i_pos]].set_non_interacting();
    } else {
      itime_vertices_[pos[i_pos]].set_interacting();
      itime_vertices_[pos[i_pos]].set_af_state(new_spins[i_pos]);
    }
  }
  sanity_check();
}

template<typename T>
void SubmatrixUpdate<T>::reject_spin_flip() {
  assert(state==TRYING_SPIN_FLIP);

  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    if (ops_rem.size()==0 && ops_ins.size()==0 && ops_replace.size()==0) {
      continue;
    } else if (ops_rem[flavor].size()==0 && ops_ins[flavor].size()>0 && ops_replace[flavor].size()==0) {
      gamma_matrices_[flavor].reject_add();
    } else if (ops_rem[flavor].size()>0 && ops_ins[flavor].size()==0 && ops_replace[flavor].size()==0) {
      gamma_matrices_[flavor].reject_remove();
    } else if (ops_rem[flavor].size()==0 && ops_ins[flavor].size()==0 && ops_replace[flavor].size()>0) {
      throw std::runtime_error("Not implemented");
      //gamma_matrices_[flavor].reject_replace();
    } else {
      throw std::runtime_error("Not implemented");
    }
  }
  sanity_check();
}
