#include "../submatrix.hpp"
/*
 * Implementation of InvAMatrix<T>
 */

template<typename T>
InvAMatrix<T>::InvAMatrix() :
    matrix_(0,0),
    creators_(0),
    annihilators_(0),
    alpha_(0),
    vertex_info_(0),
    G0_cache(0,0),
    index_G0_cache(0),
    num_entry_G0_cache(0)
{
  assert(annihilators_.size()==0);
  assert(creators_.size()==0);
}

template<typename T>
void
InvAMatrix<T>::push_back_op(const creator& cdag_op, const annihilator& c_op, T alpha, const vertex_info_type& vertex_info) {
  assert(annihilators_.size()==creators_.size());
  assert(annihilators_.size()==alpha_.size());
  assert(annihilators_.size()==vertex_info_.size());

  creators_.push_back(cdag_op);
  annihilators_.push_back(c_op);
  alpha_.push_back(alpha);
  vertex_info_.push_back(vertex_info);
}

template<typename T>
template<typename SPLINE_G0_TYPE>
void InvAMatrix<T>::extend(const SPLINE_G0_TYPE& spline_G0) {
  const int noperators = matrix_.size2();//num of operators corresponding to interacting vertices
  const int nops_add = creators_.size()-noperators;//num of operators corresponding to non-interacting vertices
  if (nops_add==0) {
    return;
  }

  //extend tilde A and set new elements to zero
  matrix_.conservative_resize(noperators + nops_add, noperators + nops_add);
  for (int j=noperators; j < noperators + nops_add; ++j) {
    for (int i=0; i < noperators + nops_add; ++i) {
      matrix_(i,j) = 0.0;
      matrix_(j,i) = 0.0;
    }
  }
  for (int i = 0; i < nops_add; ++i) {
    matrix_(i + noperators, i + noperators) = (T)1.0;
  }

  //initialize cache
  G0_cache.destructive_resize(matrix_.size1(), 100);
  index_G0_cache.resize(matrix_.size1());
  std::fill(index_G0_cache.begin(), index_G0_cache.end(), -1);
  num_entry_G0_cache = 0;

  if (noperators==0) {
    return;
  }

  //compute entries of B
  static alps::numeric::matrix<T> B;
  B.destructive_resize(nops_add, noperators);
  for (int j = 0; j < noperators; ++j) {
    for (int i = 0; i < nops_add; ++i) {
      B(i, j) = -spline_G0(annihilators_[i + noperators], creators_[j]) * (eval_f(alpha_[j]) - 1.0);
    }
  }

  //compute entries in the right lower block of A^{-1}
  //alps::numeric::submatrix_view<T> invA_view(matrix_, 0, 0, noperators, noperators);
  //alps::numeric::submatrix_view<T> B_invB_view(matrix_, noperators, 0, nops_add, noperators);
  //mygemm(-1.0, B, invA_view, (T) 0.0, B_invB_view);

  matrix_.block(noperators, 0, nops_add, noperators) = - B.block() * matrix_.block(0, 0, noperators, noperators);

  sanity_check(spline_G0);
}

template<typename T>
template<typename SPLINE_G0_TYPE>
bool InvAMatrix<T>::sanity_check(const SPLINE_G0_TYPE& spline_G0) const {
  bool result = true;
#ifndef NDEBUG
  assert(matrix_.size1()==matrix_.size2());
  assert(matrix_.size1()==creators_.size());
  assert(creators_.size()==annihilators_.size());
  assert(creators_.size()==alpha_.size());
  assert(creators_.size()==vertex_info_.size());

  result = result && (matrix_.size1()==matrix_.size2());
  result = result && (matrix_.size1()==creators_.size());
  result = result && (creators_.size()==annihilators_.size());
  result = result && (creators_.size()==alpha_.size());
  result = result && (creators_.size()==vertex_info_.size());
#endif
  return result;
}

/*
 * recompute A^{-1} and return sign(det(A)) and sign(det(1-F))
 */
template<typename T>
template<typename SPLINE_G0_TYPE>
std::pair<T,T> InvAMatrix<T>::recompute_matrix(const SPLINE_G0_TYPE& spline_G0, bool check_error) {
  const int Nv = annihilators_.size();

  if (Nv==0) return std::make_pair((T)1.0, (T)1.0);

  alps::numeric::matrix<T> matrix_bak;

  if (check_error) {
    matrix_bak = matrix_;
  }

  T sign_f_prod = 1.0;

  std::vector<T> F(Nv);
  for (int i=0; i<Nv; ++i) {
    if (alpha_at(i)==ALPHA_NON_INT) {
      throw std::logic_error("Encountered an operator corresponding to a non-interacting vertex in InvAMatrix<T>::recompute_matrix");
    }
    F[i] = eval_f(alpha_at(i));
    sign_f_prod *= (1.0-F[i])/std::abs(1.0-F[i]);
    //std::cout << "debug sign_f_prod " << i << " " << sign_f_prod << " " << F[i] << " " << alpha_at(i) << std::endl;
  }
  matrix_.conservative_resize(Nv, Nv);
  for (int j=0; j<Nv; ++j) {
    for (int i=0; i<Nv; ++i) {
      matrix_(i,j) = -spline_G0(annihilators_[i], creators_[j])*(F[j]-1.0);
    }
    matrix_(j,j) += F[j];
  }
  const T sign_det = alps::fastupdate::phase_of_determinant(matrix_);
  matrix_.invert();

  if (check_error) {
    double max_diff = -1.0, max_abs_val = 0.0;
    for (int j=0; j<Nv; ++j) {
      for (int i = 0; i < Nv; ++i) {
        max_diff = std::max(max_diff, std::abs(matrix_(i,j)-matrix_bak(i,j)));
        max_abs_val = std::max(max_abs_val, std::abs(matrix_bak(i,j)));
      }
    }
    if (max_diff/max_abs_val>1E-8) {
      std::cout << " max diff in A^{-1} is " << max_diff << ", max abs value is " << max_abs_val << " . " << std::endl;
    }
  }
  return std::make_pair(sign_det,sign_f_prod);
}

/*
 * compute det(1-F)
 */
template<typename T>
T InvAMatrix<T>::compute_f_prod() const {
  const int Nv = annihilators_.size();

  if (Nv==0) return (T)1.0;

  T f_prod = 1.0;
  for (int i=0; i<Nv; ++i) {
    f_prod *= 1.0-eval_f(alpha_at(i));
  }
  return f_prod;
}

/*
 * Implementation of InvAMatrixFlavors<T>
 */
template<typename T>
template<typename SPLINE_G0_TYPE>
void
InvAMatrixFlavors<T>::add_non_interacting_vertices(general_U_matrix<T>* p_Uijkl,
  const SPLINE_G0_TYPE& spline_G0, const itime_vertex_container& itime_vertices, int begin_index) {
  //add operators
  for (int iv=begin_index; iv<itime_vertices.size(); ++iv) {
    const itime_vertex& v = itime_vertices[iv];
    const vertex_definition<T>& vdef = p_Uijkl->get_vertex(v.type());
    assert (v.is_non_interacting());
    for (int rank=0; rank<vdef.rank(); ++rank) {
      const int flavor_rank = vdef.flavors()[rank];
      operator_time op_t(v.time(), -rank);
      sub_matrices_[flavor_rank].push_back_op(
          creator(flavor_rank, vdef.sites()[2*rank], op_t),
          annihilator(flavor_rank, vdef.sites()[2*rank+1], op_t),
          ALPHA_NON_INT,
          boost::tuple<vertex_t,size_t,my_uint64 >(v.type(), rank, v.unique_id()));
    }
  }

  //extend A^{-1}
  for (int flavor=0; flavor<sub_matrices_.size(); ++flavor) {
    sub_matrices_[flavor].extend(spline_G0);
  }
}


/*
 * Add initiliaze vectors of operators with interacting vertices. A^{-1} is then constructed.
 * Return sign(det(A)) and sign(Monte Carlo weight)=sign(det(A)*sign(-U)/det(1-F))
 */
template<typename T>
template<typename SPLINE_G0_TYPE>
std::pair<T,T>
InvAMatrixFlavors<T>::init(general_U_matrix<T>* p_Uijkl,
                                                   const SPLINE_G0_TYPE& spline_G0, const itime_vertex_container& itime_vertices, int begin_index) {
  //add operators
  T U_sign = 1.0;
  for (int iv=begin_index; iv<itime_vertices.size(); ++iv) {
    const itime_vertex& v = itime_vertices[iv];
    const vertex_definition<T>& vdef = p_Uijkl->get_vertex(v.type());
    assert (!v.is_non_interacting());
    for (int rank=0; rank<vdef.rank(); ++rank) {
      const int flavor_rank = vdef.flavors()[rank];
      operator_time op_t(v.time(), -rank);
      const T alpha = v.is_non_interacting() ? ALPHA_NON_INT : vdef.get_alpha(v.af_state(), rank);
      sub_matrices_[flavor_rank].push_back_op(
          creator(flavor_rank, vdef.sites()[2*rank], op_t),
          annihilator(flavor_rank, vdef.sites()[2*rank+1], op_t),
          alpha,
          boost::tuple<vertex_t,size_t,my_uint64 >(v.type(), rank, v.unique_id()));
    }
    U_sign *= mysign(-vdef.Uval());
  }

  T sign_det_A, f_sign;
  boost::tie(sign_det_A,f_sign) = recompute_matrix(spline_G0, false);
  T w = sign_det_A*U_sign/f_sign;
  return std::make_pair(sign_det_A, mysign(w));
}


template<typename T>
void InvAMatrix<T>::remove_rows_cols(const std::vector<int>& rows_cols) {
  const int Nv = matrix_.size1();
  const int n_rows = rows_cols.size();
  assert(n_rows<=Nv);
  for (int i=0; i<rows_cols.size(); ++i) {
    swap_rows_cols(rows_cols[n_rows-1-i], Nv-1-i);
  }
  creators_.resize(Nv-n_rows);
  annihilators_.resize(Nv-n_rows);
  alpha_.resize(Nv-n_rows);
  vertex_info_.resize(Nv-n_rows);
  matrix_.conservative_resize(Nv-n_rows, Nv-n_rows);
}

template<typename T>
void
InvAMatrix<T>::remove_rows_cols(const std::vector<my_uint64>& v_uid) {
  std::vector<int> rows_removed, tmp;
  for (int it=0; it<v_uid.size(); ++it) {
    tmp = find_row_col(v_uid[it]);
    for (std::vector<int>::iterator it=tmp.begin(); it!=tmp.end(); ++it)
      rows_removed.push_back(*it);
  }
  if (rows_removed.size()==0) return;
  std::sort(rows_removed.begin(), rows_removed.end());
  remove_rows_cols(rows_removed);
}

template<typename T>
template<typename SPLINE_G0_TYPE>
void
InvAMatrix<T>::update_matrix(const InvGammaMatrix<T>& inv_gamma, const SPLINE_G0_TYPE& spline_G0) {

  const int nop = inv_gamma.matrix().size1();
  const int N = matrix_.size2();

  sanity_check(spline_G0);

  if (nop==0) return;

  invA0.destructive_resize(nop, N);
  G0_left.destructive_resize(N, nop);
  G0_inv_gamma.destructive_resize(N, nop);

  pl.resize(nop);
  for (int l=0; l<nop; ++l) {
    pl[l] = inv_gamma.pos_in_invA(l);
  }
  for (int l=0; l<nop; ++l) {
    invA0.row(l) = matrix_.row(pl[l]);
  }
  for (int l=0; l<nop; ++l) {
    auto view = G0_left.block(0, l, N, 1);
    eval_Gij_col(spline_G0, pl[l], view);
  }

  G0_inv_gamma.block() = - G0_left.block() * inv_gamma.matrix().block();
  matrix_.block() += G0_inv_gamma.block() * invA0.block();

  for (int l=0; l < nop; ++l) {
    assert(inv_gamma.alpha(l)!=inv_gamma.alpha0(l));
    if (inv_gamma.alpha(l)==inv_gamma.alpha0(l)) {
      throw std::runtime_error("fatal error in update_matrix.");
    }

    const int i_row = pl[l];
    assert(inv_gamma.alpha0(l)==alpha_at(i_row));
    const T coeff = 1.0/(
        1.0 +
        gamma_func(
            eval_f(inv_gamma.alpha(l)),
            eval_f(inv_gamma.alpha0(l))
        )
    );
    for (int i_col=0; i_col<N; ++i_col) {
      matrix_(i_row, i_col) *= coeff;
    }
  }

  //update alpha
  for (int l=0; l<nop; ++l) {
    assert(pl[l]>=0 && pl[l]<N);
    alpha_[pl[l]] = inv_gamma.alpha(l);
  }

  //clear cache
  G0_cache.destructive_resize(0,0);
  index_G0_cache.resize(0);
  num_entry_G0_cache = 0;
}

//compute M=(G-alpha)^-1 from A^-1
// using the relation M = (1-f) A^-1
template<typename T>
template<typename SPLINE_G0_TYPE>
void InvAMatrix<T>::compute_M(alps::numeric::matrix<T>& M, const SPLINE_G0_TYPE& spline_G0) const {
  const int N = matrix_.size2();
  M.destructive_resize(N, N);

  if (N==0) return;

  std::vector<T> coeff(N);
  for (int i=0; i<N; ++i) {
    coeff[i] = 1.0-eval_f(alpha_at(i));
  }

  for (int j=0; j<N; ++j) {
    for (int i=0; i<N; ++i) {
      M(i,j) = coeff[i]*matrix_(i,j);
    }
  }

#ifndef NDEBUG
  alps::numeric::matrix<T> G0(N,N), G0_M(N,N);
  for (int j=0; j<N; ++j) {
    for (int i=0; i<N; ++i) {
      G0(i,j) = spline_G0(annihilators()[i], creators()[j]);
    }
    G0(j,j) -= alpha_at(j);
  }
  //mygemm((T) 1.0, G0, M, (T) 0.0, G0_M);
  G0_M.block() = G0.block() * M.block();
  for (int j=0; j<N; ++j) {
    for (int i = 0; i < N; ++i) {
      if (i==j) {
        assert(std::abs(G0_M(i,j)-1.0)<1E-5);
      } else {
        assert(std::abs(G0_M(i,j))<1E-5);
      }
    }
  }
#endif
}

// G_{ij} = sum_p (A^{-1})_{ip}, G0_{pj}
// cols specifies {j}
template<typename T>
template<typename SPLINE_G0_TYPE, typename M>
void InvAMatrix<T>::eval_Gij_col(const SPLINE_G0_TYPE& spline_G0, int col, M& Gij) const {
  //static alps::numeric::matrix<T> G0;
 assert (col>=0);

  const int Nv = matrix_.size1();
  const T alpha_col = alpha_at(col);
  //assert(Gij.size1()==Nv);
  //assert(Gij.size2()==1);
  if (alpha_col!=ALPHA_NON_INT) {
    const T fj = eval_f(alpha_col);
    for (int iv=0; iv<Nv; ++iv) {
      Gij(iv,0) = (fj*matrix_(iv,col))/(fj-1.0);
    }
    Gij(col,0) = (fj*matrix_(col,col)-1.0)/(fj-1.0);
  } else {
    auto G0_view = compute_G0_col(spline_G0, col);
    //mygemm((T) 1.0, matrix_, G0_view, (T) 0.0, Gij);
    Gij = matrix_.block() * G0_view;
  }
}

// Compute G_{ij} = sum_p (A^{-1})_{ip}, G0_{pj} for given sets of i and j.
template<typename T>
template<typename SPLINE_G0_TYPE, typename M>
void InvAMatrix<T>::eval_Gij_cols_rows(const SPLINE_G0_TYPE& spline_G0,
                                       const std::vector<int>& rows,
                                       const std::vector<int>& cols, M& result) const {
  const int Nv = matrix_.size1();
  const int n_rows = rows.size();

  std::vector<int> idx_cols_slow, cols_slow;
  for (int icol = 0; icol < cols.size(); ++icol) {
    int col = cols[icol];
    T alpha_col = alpha_at(col);
    if (alpha_col!=ALPHA_NON_INT) {
      // Use Fast formula
      // Eq. (A3) in Nomura et al (2014)
      auto view = result.block(0, icol, n_rows, 1);
      eval_Gij_col_part(spline_G0, rows, col, view);

      /*
      for (int iv=0; iv<n_rows; ++iv) {
        std::cout << " iv " << iv << " " << result(iv,icol) << std::endl;
      }

      const T fj = eval_f(alpha_col);
      for (int iv=0; iv<n_rows; ++iv) {
        if (rows[iv] != col) {
          //result(iv,icol) = (fj*matrix_(rows[iv],col))/(fj-1.0);
          auto r = (fj*matrix_(rows[iv],col))/(fj-1.0);
          std::cout << " iv2 " << iv << " " << r << std::endl;
        } else {
          auto r = (fj*matrix_(rows[iv],col)-1.0)/(fj-1.0);
          std::cout << " iv2 " << iv << " " << r << std::endl;
          //result(iv,icol) = (fj*matrix_(rows[iv],col)-1.0)/(fj-1.0);
        }
      }
      */

      //for (int iv=0; iv<n_rows; ++iv) {
        //std::cout << " iv2 " << iv << " " << result(iv,icol) << std::endl;
      //}

    } else {
      // Use Slow formula later
      idx_cols_slow.push_back(icol);
      cols_slow.push_back(col);
    }
  }

  // Use Slow Formula with optimization of blas3 level
  // Eq. (A4) in Nomura et al (2014)
  using matrix_type =  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
  matrix_type G0_tmp(Nv, cols_slow.size());
  for (int icol = 0; icol < cols_slow.size(); ++icol) {
    G0_tmp.block(0, icol, Nv, 1) = compute_G0_col(spline_G0, cols_slow[icol]);
  }

  matrix_type invA_tmp(rows.size(), Nv);
  for (int iv=0; iv<n_rows; ++iv) {
    invA_tmp.block(iv, 0, 1, Nv) = matrix_.block(rows[iv], 0, 1, Nv);
  }

  matrix_type Gij_slow = (invA_tmp * G0_tmp);

  for (int icol = 0; icol < cols_slow.size(); ++icol) {
    result.block(0, idx_cols_slow[icol], n_rows, 1) = Gij_slow.block(0, icol, n_rows, 1);
  }

  //for (int iv=0; iv<n_rows; ++iv) {
    //Gij(iv,0) = (matrix_.block(rows[iv], 0, 1, Nv) * G0_view)(0,0);
  //}

};

// G_{ij} = sum_p (A^{-1})_{ip}, G0_{pj}
// cols specifies {j}
template<typename T>
template<typename SPLINE_G0_TYPE, typename M>
void InvAMatrix<T>::eval_Gij_col_part(const SPLINE_G0_TYPE& spline_G0, const std::vector<int>& rows, int col, M& Gij) const {
  assert (col>=0);

  const int Nv = matrix_.size1();
  const T alpha_col = alpha_at(col);
  const int n_rows = rows.size();

  assert(Gij.rows()==n_rows);
  assert(Gij.cols()==1);

  if (alpha_col!=ALPHA_NON_INT) {
    const T fj = eval_f(alpha_col);
    for (int iv=0; iv<n_rows; ++iv) {
      if (rows[iv]!=col) {
        Gij(iv,0) = (fj*matrix_(rows[iv],col))/(fj-1.0);
      } else {
        Gij(iv,0) = (fj*matrix_(rows[iv],col)-1.0)/(fj-1.0);
      }
    }
  } else {
    throw std::runtime_error("Do not use this anymore");
    alps::numeric::submatrix_view<T> G0_view = compute_G0_col(spline_G0, col);
    alps::numeric::matrix<T> G_tmp(1,1);
    for (int iv=0; iv<n_rows; ++iv) {
      //alps::numeric::submatrix_view<T> invA_view(matrix_, rows[iv], 0, 1, Nv);
      //auto invA_view = matrix_.block(rows[iv], 0, 1, Nv);
      //mygemm((T) 1.0, invA_view, G0_view, (T) 0.0, G_tmp);
      //G_tmp.block() = matrix_.block(rows[iv], 0, 1, Nv) * G0_view;
      Gij(iv,0) = (matrix_.block(rows[iv], 0, 1, Nv) * G0_view)(0,0);
    }

  }
}

template<typename T>
template<typename SPLINE_G0_TYPE>
alps::numeric::submatrix_view<T> InvAMatrix<T>::compute_G0_col(const SPLINE_G0_TYPE& spline_G0, int col) const {
  assert(col>=0);
  //look up cache
  const int Nv = matrix_.size1();
  assert(col<index_G0_cache.size());
  if (index_G0_cache[col]>=0) {
    assert(index_G0_cache[col]<G0_cache.size2());
    return G0_cache.block(0, index_G0_cache[col], Nv, 1);
  } else {
    const int index = num_entry_G0_cache;
    index_G0_cache[col] = num_entry_G0_cache;
    if (G0_cache.size2()<=num_entry_G0_cache) {
      G0_cache.conservative_resize(Nv, static_cast<int>(1.5*num_entry_G0_cache)+1);
    }
    for (int iv=0; iv<Nv; ++iv) {
      G0_cache(iv,index) = spline_G0(annihilators_[iv], creators_[col]);
    }
    ++num_entry_G0_cache;
    return G0_cache.block(0, index, Nv, 1);
  }
};

/*
 * Implementation of InvAMatrixFlavors<T>
 */
template<typename T>
template<typename SPLINE_G0_TYPE>
void
InvAMatrixFlavors<T>::update_matrix(const std::vector<InvGammaMatrix<T> >& inv_gamma_flavors, const SPLINE_G0_TYPE& spline_G0) {
  for (int flavor=0; flavor<sub_matrices_.size(); ++flavor) {
    sub_matrices_[flavor].update_matrix(inv_gamma_flavors[flavor], spline_G0);
  }
}

template<typename T>
void
InvAMatrixFlavors<T>::remove_rows_cols(const std::vector<my_uint64>& v_uid) {
  for (int flavor=0; flavor<sub_matrices_.size(); ++flavor) {
    sub_matrices_[flavor].remove_rows_cols(v_uid);
  }
}

template<typename T>
template<typename SPLINE_G0_TYPE>
bool
InvAMatrixFlavors<T>::sanity_check(const SPLINE_G0_TYPE& spline_G0, general_U_matrix<T>* p_Uijkl,
                                   const itime_vertex_container& itime_vertices) const {
  bool result = true;
#ifndef NDEBUG
  const int n_flavors = sub_matrices_.size();
  std::vector<std::vector<creator> > creators_scr(n_flavors);
  std::vector<std::vector<annihilator> > annihilators_scr(n_flavors);
  std::vector<std::vector<T> > alpha_scr(n_flavors);
  std::vector<std::vector<my_uint64 > > v_uid_scr(n_flavors);
  std::vector<std::vector<my_uint64 > > rank_scr(n_flavors);

  for (int iv=0; iv<itime_vertices.size(); ++iv) {
    const itime_vertex& v = itime_vertices[iv];
    const vertex_definition<T>& vdef = p_Uijkl->get_vertex(v.type());
    for (int rank=0; rank<v.rank(); ++rank) {
      const int flavor_rank = vdef.flavors()[rank];
      operator_time op_t(v.time(), -rank);
      creators_scr[flavor_rank].push_back(
          creator(flavor_rank, vdef.sites()[2*rank], op_t)
      );
      annihilators_scr[flavor_rank].push_back(
          annihilator(flavor_rank, vdef.sites()[2*rank+1], op_t)
      );
      alpha_scr[flavor_rank].push_back(vdef.get_alpha(v.af_state(), rank));
      v_uid_scr[flavor_rank].push_back(v.unique_id());
      rank_scr[flavor_rank].push_back(rank);
    }
  }

  for (int flavor=0; flavor<n_flavors; ++flavor) {
    result = result && sub_matrices_[flavor].sanity_check(spline_G0);
    assert(sub_matrices_[flavor].annihilators().size()==annihilators_scr[flavor].size());
    assert(sub_matrices_[flavor].creators().size()==creators_scr[flavor].size());

    //check operators one by one
    for (int iop=0; iop<creators_scr[flavor].size(); ++iop) {
      //const int pos = sub_matrices_[flavor].find_row_col(creators_scr[flavor][iop].t());
      const int pos = sub_matrices_[flavor].find_row_col(v_uid_scr[flavor][iop], rank_scr[flavor][iop]);
      assert(sub_matrices_[flavor].creators()[pos]==creators_scr[flavor][iop]);
      assert(sub_matrices_[flavor].annihilators()[pos]==annihilators_scr[flavor][iop]);
      assert(sub_matrices_[flavor].alpha_at(pos)==alpha_scr[flavor][iop]);
    }
  }
#endif
  return result;
}

/*
 * recompute A^{-1} and return sign(det(A)) and sign(det(1-F))
 */
template<typename T>
template<typename SPLINE_G0_TYPE>
std::pair<T,T> InvAMatrixFlavors<T>::recompute_matrix(const SPLINE_G0_TYPE& spline_G0, bool check_error) {
  T sign_detA = 1.0, f_sign = 1.0;
  for (int flavor=0; flavor<sub_matrices_.size(); ++flavor) {
    T sign_detA_tmp, f_sign_tmp;
    boost::tie(sign_detA_tmp,f_sign_tmp) = sub_matrices_[flavor].recompute_matrix(spline_G0, check_error);
    sign_detA *= sign_detA_tmp;
    f_sign *= f_sign_tmp;
  }
  return std::make_pair(sign_detA,f_sign);
}

