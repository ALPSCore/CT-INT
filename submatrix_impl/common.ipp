#include "../submatrix.hpp"

template<typename T>
SubmatrixUpdate<T>::SubmatrixUpdate(int k_ins_max, int n_flavors, SPLINE_G0_TYPE spline_G0, general_U_matrix<T>* p_Uijkl, double beta) :
    k_ins_max_(k_ins_max),
    spline_G0_(spline_G0),
    p_Uijkl_(p_Uijkl),
    beta_(beta),
    coeff_det(n_flavors%2==0 ? 1 : -1),
    alpha_scale_(1),
    invA_(n_flavors),
    sign_det_A_(1.0),
    sign_(1.0),
    state(READY_FOR_UPDATE),
    gamma_matrices_(n_flavors),
    ops_rem(n_flavors),
    ops_ins(n_flavors),
    ops_replace(n_flavors),
    current_vertex_id_(0)
{
  for (int i=0; i<101; ++i) {
    std::cout << " i " << i << " " <<
        spline_G0(
            annihilator(0,0,operator_time(beta*i/100.0-beta, 0)),
            creator(0,0,operator_time(0.0, 0))) << " " <<
        spline_G0(
        annihilator(0,0,operator_time(beta*i/100.0, 0)),
        creator(0,0,operator_time(0.0, 0)))<< " " <<
    spline_G0(
        annihilator(0,0,operator_time(beta*i/100.0+beta, 0)),
        creator(0,0,operator_time(0.0, 0))

    ) << std::endl;
    //std::cout << " i " << 0 << " " <<spline_G0(
        //annihilator(0,0,operator_time(0.0, 0)),
        //creator(0,0,operator_time(0.0, 0))) << std::endl;
  }
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
    sign_det_A_(1.0),
    sign_(1.0),
    state(READY_FOR_UPDATE),
    gamma_matrices_(n_flavors),
    ops_rem(n_flavors),
    ops_ins(n_flavors),
    ops_replace(n_flavors),
    current_vertex_id_(0)
{
  if (itime_vertices_init.size()!=0) {

    itime_vertices_ = itime_vertices_init;
    for (int iv=0; iv<itime_vertices_.size(); ++iv) {
      itime_vertices_[iv].set_unique_id(gen_new_vertex_id());
      //std::cout << "uid " << itime_vertices_[iv].unique_id() << std::endl;
    }
    boost::tie(sign_det_A_,sign_) = invA_.init(p_Uijkl, spline_G0, itime_vertices_, 0);


    //sign_ = 1.0;
    //for (int iv=0; iv<itime_vertices_init.size(); ++iv) {
      //sign_ *= mysign(p_Uijkl->get_vertex(itime_vertices_init[iv].type()).Uval());
    //}
    //std::vector<alps::numeric::matrix<T> > M(n_flavors);
    //compute_M_from_scratch(M);
    //for (int flavor=0; flavor<n_flavors; ++flavor) {
      //const T det_M = alps::numeric::safe_determinant(M[flavor]);
      //sign_ *= det_M/std::abs(det_M);
    //}

    std::cout << "Reconstructed sign_det_A, sign " << sign_det_A_ << " " << sign_ << std::endl;
  }
  /*
  for (int i=0; i<101; ++i) {
    std::cout << " i " << i << " " <<
    spline_G0(
        annihilator(0,0,operator_time(beta*i/100.0-beta, 0)),
        creator(0,0,operator_time(0.0, 0))) << " " <<
        spline_G0(
            annihilator(0,0,operator_time(beta*i/100.0, 0)),
            creator(0,0,operator_time(0.0, 0))) << " " <<
            spline_G0(
                annihilator(0,0,operator_time(beta*i/100.0+beta, 0)),
                creator(0,0,operator_time(0.0, 0))) << std::endl;
  }
  */
}

template<typename T, typename SPLINE_G0_TYPE>
T eval_Gij(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0, int row_A, int col_A) {
  static alps::numeric::matrix<T> G0, invA_G0(1,1);

  const T alpha_col = invA.alpha_at(col_A);
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

/*
// G_{ij} = sum_p (A^{-1})_{ip}, G0_{pj}
// cols specifies {j}
template<typename T, typename SPLINE_G0_TYPE, typename M>
void eval_Gij_col(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0, int col, M& Gij) {
  static alps::numeric::matrix<T> G0;

  const int Nv = invA.matrix().num_rows();
  const T alpha_col = invA.alpha_at(col);
  assert(Gij.num_rows()==Nv);
  assert(Gij.num_cols()==1);
  if (alpha_col!=ALPHA_NON_INT) {
    const T fj = eval_f(alpha_col);
    for (int iv=0; iv<Nv; ++iv) {
      Gij(iv,0) = (fj*invA.matrix()(iv,col))/(fj-1.0);
    }
    Gij(col,0) = (fj*invA.matrix()(col,col)-1.0)/(fj-1.0);
  } else {
    G0.resize_values_not_retained(Nv,1);
    for (int iv=0; iv<Nv; ++iv) {
      G0(iv,0) = spline_G0(invA.annihilators()[iv], invA.creators()[col]);
    }
    mygemm((T) 1.0, invA.matrix(), G0, (T) 0.0, Gij);
  }
}

// G_{ij} = sum_p (A^{-1})_{ip}, G0_{pj}
// cols specifies {j}
template<typename T, typename SPLINE_G0_TYPE, typename M>
void eval_Gij_col(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0, const std::vector<int>& rows, int col, M& Gij) {
  static alps::numeric::matrix<T> G0, invA_tmp;

  const int Nv = invA.matrix().num_rows();
  const T alpha_col = invA.alpha_at(col);
  const int n_rows = rows.size();

  assert(Gij.num_rows()==n_rows);
  assert(Gij.num_cols()==1);

  if (alpha_col!=ALPHA_NON_INT) {
    const T fj = eval_f(alpha_col);
    for (int iv=0; iv<n_rows; ++iv) {
      if (rows[iv]!=col) {
        Gij(iv,0) = (fj*invA.matrix()(rows[iv],col))/(fj-1.0);
      } else {
        Gij(iv,0) = (fj*invA.matrix()(rows[iv],col)-1.0)/(fj-1.0);
      }
    }
  } else {
    G0.resize_values_not_retained(Nv,1);
    invA_tmp.resize_values_not_retained(n_rows,Nv);
    for (int iv=0; iv<n_rows; ++iv) {
      for (int iv2=0; iv2<Nv; ++iv2) {
        invA_tmp(iv,iv2) = invA.matrix()(rows[iv], iv2);
      }
    }
    for (int iv=0; iv<Nv; ++iv) {
      G0(iv,0) = spline_G0(invA.annihilators()[iv], invA.creators()[col]);
    }
    mygemm((T) 1.0, invA_tmp, G0, (T) 0.0, Gij);
  }
}
*/

template<typename T>
bool SubmatrixUpdate<T>::sanity_check() {
  bool result = true;
#ifndef NDEBUG
  //check gamma^{-1}
  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    result = result && gamma_matrices_[flavor].sanity_check(invA_[flavor], spline_G0_);
  }

  //check A^{-1}
  if (state==READY_FOR_UPDATE) {
    result = result && invA_.sanity_check(spline_G0_, p_Uijkl_, itime_vertices_);

    const T sign_det_A_bak = sign_det_A_;
    const T sign_bak = sign_;
    recompute_matrix(true);//recompute A^{-1}
    assert(my_equal(sign_bak, sign_, 1E-5));
    assert(my_equal(sign_det_A_bak, sign_det_A_, 1E-5));
  }
#endif
  return result;
}

/*
 * Recompute A^{-1} and sign of Monte Carl weight.
 */
template<typename T>
void SubmatrixUpdate<T>::recompute_matrix(bool check_error) {
  if (state==READY_FOR_UPDATE) {
    const T sign_det_A_bak = sign_det_A_;
    const T sign_bak = sign_;
    T f_sign;
    boost::tie(sign_det_A_,f_sign) = invA_.recompute_matrix(spline_G0_, check_error);
    sign_ = sign_det_A_/f_sign;
    for (int iv=0; iv<itime_vertices_.size(); ++iv) {
      assert (!itime_vertices_[iv].is_non_interacting());
      sign_ *= mysign(-p_Uijkl_->get_vertex(itime_vertices_[iv].type()).Uval());
    }

    if (check_error) {
      if (!my_equal(sign_,sign_bak)) {
        std::cout << " Error in sign is " << std::abs(sign_-sign_bak) << std::endl;
        std::cout << "sign_ " << sign_ << std::endl;
        std::cout << "sign_bak " << sign_bak << std::endl;
        std::cout << "sign_det_A_ " << sign_det_A_ << std::endl;
        std::cout << "f_sign " << f_sign << std::endl;
      }
      if (!my_equal(sign_det_A_,sign_det_A_bak)) {
        std::cout << " Error in sign_det_A is " << std::abs(sign_det_A_-sign_det_A_bak) << std::endl;
      }
    }
  } else {
    throw std::logic_error("The call to SubmatrixUpdate<T>::recompute_matrix is illegal.");
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
  std::vector<my_uint64> uid_removed;
  itime_vertex_container new_itime_vertices;
  for (int iv=0; iv<itime_vertices_.size(); ++iv) {
    if (itime_vertices_[iv].is_non_interacting()) {
      uid_removed.push_back(itime_vertices_[iv].unique_id());
    } else {
      new_itime_vertices.push_back(itime_vertices_[iv]);
    }
  }
  invA_.remove_rows_cols(uid_removed);
  std::swap(itime_vertices_, new_itime_vertices);

  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    gamma_matrices_[flavor].clear();
  }

  //some check
  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    const int Nv = invA_[flavor].annihilators().size();
    if (invA_[flavor].creators().size()!=Nv) {
      throw std::logic_error("creators().size() != annihilators().size()");
    }
    for (int iv=0; iv<Nv; ++iv) {
      if (invA_[flavor].alpha_at(iv)==ALPHA_NON_INT) {
        throw std::logic_error("Found an operator corresponding to a non-interacting vertex.");
      }
    }
  }
  state = READY_FOR_UPDATE;
}

//returns the determinant ratios of A and 1-f, respectively.
template<typename T>
boost::tuple<T,T,T>
SubmatrixUpdate<T>::try_spin_flip(const std::vector<int>& pos, const std::vector<int>& new_spins) {
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

  T f_rat = 1.0, U_rat = 1.0;

  //figure out which rows/cols are to be updated
  for (int ipos=0; ipos<pos.size(); ++ipos) {
    const int iv = pos[ipos];
    const int new_spin = new_spins[ipos];
    assert(iv<itime_vertices0_.size());
    //const itime_vertex& v = itime_vertices0_[iv];
    const itime_vertex& v = itime_vertices_[iv];
    const vertex_definition<T>& vdef = p_Uijkl_->get_vertex(v.type());

    if(v.is_non_interacting() && new_spin!=NON_INT_SPIN_STATE) {
      U_rat *= -vdef.Uval();
    } else if(!v.is_non_interacting() && new_spin==NON_INT_SPIN_STATE) {
      U_rat /= -vdef.Uval();
    }

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
      const int pos_in_A = invA_[flavor_rank].find_row_col(v.unique_id(), rank);
      assert(pos_in_A>=0);
      OperatorToBeUpdated<T> elem(op_t, pos_in_A, alpha0, alpha_current, alpha_new);

      const bool already_in_gamma = (alpha0 != alpha_current);
      const bool stay_in_gamma = (alpha0 != alpha_new);

      if (already_in_gamma && !stay_in_gamma) {
        ops_rem[flavor_rank].push_back(elem);
      } else if (!already_in_gamma && stay_in_gamma) {
        ops_ins[flavor_rank].push_back(elem);
      } else if (already_in_gamma && stay_in_gamma) {
        if (alpha_current != alpha_new) {
          ops_replace[flavor_rank].push_back(elem);
        }
      } else {
        //nothing to do with Gamma
      }

      if (alpha_current!=alpha_new) {
        if (alpha_new!=ALPHA_NON_INT) {
          f_rat *= 1.0-eval_f(alpha_new);
        }
        if (alpha_current!=ALPHA_NON_INT) {
          f_rat /= 1.0-eval_f(alpha_current);
        }
      }
    }
  }

  //compute determinant ratio of A (not M):
  det_rat_A = 1.0;
  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    if (ops_rem.size()==0 && ops_ins.size()==0 && ops_replace.size()==0) {
      continue;
    } else if (ops_rem[flavor].size()==0 && ops_ins[flavor].size()>0 && ops_replace[flavor].size()==0) {
        det_rat_A *= gamma_matrices_[flavor].try_add(invA_[flavor], spline_G0_, ops_ins[flavor]);
    } else if (ops_rem[flavor].size()>0 && ops_ins[flavor].size()==0 && ops_replace[flavor].size()==0) {
      det_rat_A *= gamma_matrices_[flavor].try_remove(invA_[flavor], spline_G0_, ops_rem[flavor]);
    } else if (ops_rem[flavor].size()>0 && ops_ins[flavor].size()>0 && ops_replace[flavor].size()==0) {
      det_rat_A *= gamma_matrices_[flavor].try_add_remove(invA_[flavor], spline_G0_, ops_ins[flavor], ops_rem[flavor]);
    } else if (ops_rem[flavor].size()==0 && ops_ins[flavor].size()==0 && ops_replace[flavor].size()>0) {
      throw std::runtime_error("Not implemented: try_replace");
      //det_rat_A *= gamma_matrices_[flavor].try_replace(p_Uijkl_, invA_[flavor], ops_ins[flavor]);
    } else {
      std::cout << "debug info " << ops_rem[flavor].size() << " " << ops_ins[flavor].size() << " " << ops_replace[flavor].size() << std::endl;
      throw std::runtime_error("Not implemented: try_***");
    }
  }

  sign_rat = mysign(det_rat_A)*mysign(U_rat)/mysign(f_rat);

  return boost::make_tuple(det_rat_A, 1.0/f_rat, U_rat);
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
    } else if (ops_rem[flavor].size()>0 && ops_ins[flavor].size()>0 && ops_replace[flavor].size()==0) {
      gamma_matrices_[flavor].perform_add_remove();
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

  //std::cout << "old det_A = " << det_A_ << std::endl;

  sign_det_A_ *= mysign(det_rat_A);
  sign_ *= sign_rat;

  //if (mycast<double>(sign_)<0.0) throw std::runtime_error("Sign is negative");

  //std::cout << "using det_rat_A = " << det_rat_A << std::endl;
  //std::cout << "new det_A = " << det_A_ << std::endl;
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
    } else if (ops_rem[flavor].size()>0 && ops_ins[flavor].size()>0 && ops_replace[flavor].size()==0) {
      gamma_matrices_[flavor].reject_add_remove();
    } else if (ops_rem[flavor].size()==0 && ops_ins[flavor].size()==0 && ops_replace[flavor].size()>0) {
      throw std::runtime_error("Not implemented");
      //gamma_matrices_[flavor].reject_replace();
    } else {
      throw std::runtime_error("Not implemented");
    }
  }
  sanity_check();
}

template<typename  T>
void SubmatrixUpdate<T>::compute_M(std::vector<alps::numeric::matrix<T> >& M) {
  assert(M.size()==n_flavors());

  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    invA_[flavor].compute_M(M[flavor], spline_G0_);
  }
}

template<typename  T>
void SubmatrixUpdate<T>::compute_M_from_scratch(std::vector<alps::numeric::matrix<T> >& M) {
  assert(M.size()==n_flavors());

  M.resize(n_flavors());
  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    const int Nv = invA_[flavor].matrix().num_cols();
    M[flavor].resize(Nv, Nv);
    if (Nv>0) {
      for (int j=0; j<Nv; ++j) {
        for (int i=0; i<Nv; ++i) {
          M[flavor](i,j) = spline_G0_(invA_[flavor].annihilators()[i], invA_[flavor].creators()[j]);
        }
        M[flavor](j,j) -= invA_[flavor].alpha_at(j);
      }
      alps::numeric::inverse_in_place(M[flavor]);
    }
  }
}


/*
template<typename  T>
T SubmatrixUpdate<T>::recompute_sign(bool check_error) {
  assert(state==READY_FOR_UPDATE);
  // -U |G0-alpha| = -U |A|/|1-f|

  T prod_f = 1.0;
  for (int flavor=0; flavor<n_flavors(); ++flavor) {
    prod_f *= invA_[flavor].compute_f_prod();
  }

  T prod_U = 1.0;
  for (int iv=0; iv<itime_vertices_.size(); ++iv) {
    assert (!itime_vertices_[iv].is_non_interacting());
    prod_U *= -p_Uijkl_->get_vertex(itime_vertices_[iv].type()).Uval();
  }

  const T weight = det_A_*prod_U/prod_f;
  //std::cout << "debug_weight " << det_A_ << " " << prod_U << " " << 1.0/prod_f << std::endl;
#ifndef NDEBUG
  if (check_error) {
    //std::cout << "debug sign " << sign_<< " " << weight/std::abs(weight) << std::endl;
    assert(std::abs(sign_ - weight/std::abs(weight))<1E-5);
  }
#endif
  sign_ = weight/std::abs(weight);
}
*/

template<typename  T>
my_uint64 SubmatrixUpdate<T>::gen_new_vertex_id() {
  ++current_vertex_id_;
  return current_vertex_id_;
}
