#include "../submatrix.hpp"
/*
 * Implementation of InvAMatrix<T>
 */

template<typename T>
InvAMatrix<T>::InvAMatrix() : alpha_scale_(1.0), matrix_(0,0), creators_(0), annihilators_(0), alpha_(0), vertex_info_(0) {
  assert(annihilators_.size()==0);
  assert(creators_.size()==0);
}

template<typename T>
int
InvAMatrix<T>::find_row_col(double time) const {
  for(std::size_t i=0; i<creators_.size(); ++i) {
    if (time==creators_[i].t().time()) {
      assert(annihilators_[i].t().time()==time);
      return i;
    }
  }
  return -1;
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
  const int noperators = matrix_.num_cols();//num of operators corresponding to interacting vertices
  const int nops_add = creators_.size()-noperators;//num of operators corresponding to non-interacting vertices
  if (nops_add==0) {
    return;
  }

  //extend tilde A
  matrix_.resize(noperators + nops_add, noperators + nops_add, (T)0.0);
  for (int i = 0; i < nops_add; ++i) {
    matrix_(i + noperators, i + noperators) = (T)1.0;
  }

  if (noperators==0) {
    return;
  }

  //compute entries of B
  static alps::numeric::matrix<T> B;
  B.resize_values_not_retained(nops_add, noperators);
  for (int j = 0; j < noperators; ++j) {
    for (int i = 0; i < nops_add; ++i) {
      B(i, j) = -spline_G0(annihilators_[i + noperators], creators_[j]) * (eval_f(alpha_[j]) - 1.0);
    }
  }

  //compute entries in the right lower block of A^{-1}
  alps::numeric::submatrix_view<T> invA_view(matrix_, 0, 0, noperators, noperators);
  alps::numeric::submatrix_view<T> B_invB_view(matrix_, noperators, 0, nops_add, noperators);
  mygemm(-1.0, B, invA_view, (T) 0.0, B_invB_view);

  sanity_check(spline_G0);
}

template<typename T>
template<typename SPLINE_G0_TYPE>
void InvAMatrix<T>::sanity_check(const SPLINE_G0_TYPE& spline_G0) const {
#ifndef NDEBUG
  assert(num_rows(matrix_)==num_cols(matrix_));
  assert(num_rows(matrix_)==creators_.size());
  assert(creators_.size()==annihilators_.size());
  assert(creators_.size()==alpha_.size());
  assert(creators_.size()==vertex_info_.size());
#endif
}

template<typename T>
template<typename SPLINE_G0_TYPE>
void InvAMatrix<T>::recompute(const SPLINE_G0_TYPE& spline_G0, bool check_error) {
  const int Nv = annihilators_.size();

  if (Nv==0) return;

  alps::numeric::matrix<T> matrix_bak;

  if (check_error) {
    matrix_bak = matrix_;
    assert(matrix_.num_cols()==Nv);
    assert(matrix_.num_rows()==Nv);
  }

  std::vector<T> F(Nv);
  for (int i=0; i<Nv; ++i) {
    F[i] = eval_f(alpha_at(i));
  }
  matrix_.resize(Nv, Nv);
  for (int j=0; j<Nv; ++j) {
    for (int i=0; i<Nv; ++i) {
      matrix_(i,j) = -spline_G0(annihilators_[i], creators_[j])*(F[j]-1.0);
    }
    matrix_(j,j) += F[j];
  }
  alps::numeric::inverse_in_place(matrix_);

  if (check_error) {
    double max_diff = -1.0, max_abs_val = 0.0;
    for (int j=0; j<Nv; ++j) {
      for (int i = 0; i < Nv; ++i) {
        max_diff = std::max(max_diff, std::abs(matrix_(i,j)-matrix_bak(i,j)));
        max_abs_val = std::max(max_abs_val, std::abs(matrix_bak(i,j)));
      }
    }
    if (max_diff>1E-8) {
      std::cout << " max diff in A^{-1} is " << max_diff << ", max abs value is " << max_abs_val << " . " << std::endl;
    }
  }
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
          std::pair<vertex_t,size_t>(v.type(), rank));
    }
  }

  //extend A^{-1}
  for (int flavor=0; flavor<sub_matrices_.size(); ++flavor) {
    sub_matrices_[flavor].extend(spline_G0);
  }
}


template<typename T>
template<typename SPLINE_G0_TYPE>
void
InvAMatrixFlavors<T>::add_interacting_vertices(general_U_matrix<T>* p_Uijkl,
                                                   const SPLINE_G0_TYPE& spline_G0, const itime_vertex_container& itime_vertices, int begin_index) {
  //add operators
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
          std::pair<vertex_t,size_t>(v.type(), rank));
    }
  }

  for (int flavor=0; flavor<sub_matrices_.size(); ++flavor) {
    sub_matrices_[flavor].recompute(spline_G0, false);
  }
}


template<typename T>
void InvAMatrix<T>::remove_rows_cols(const std::vector<int>& rows_cols) {
  const int Nv = matrix_.num_rows();
  const int n_rows = rows_cols.size();
  assert(n_rows<=Nv);
  for (int i=0; i<rows_cols.size(); ++i) {
    swap_rows_cols(rows_cols[n_rows-1-i], Nv-1-i);
  }
  creators_.resize(Nv-n_rows);
  annihilators_.resize(Nv-n_rows);
  alpha_.resize(Nv-n_rows);
  vertex_info_.resize(Nv-n_rows);
  matrix_.resize(Nv-n_rows, Nv-n_rows);
}

template<typename T>
void
InvAMatrix<T>::remove_rows_cols(const std::vector<double>& times) {
  std::vector<int> rows_removed;
  for (int it=0; it<times.size(); ++it) {
    rows_removed.push_back(find_row_col(times[it]));
  }
  if (rows_removed.size()==0) return;
  std::sort(rows_removed.begin(), rows_removed.end());
  remove_rows_cols(rows_removed);
}

template<typename T>
template<typename SPLINE_G0_TYPE>
void
InvAMatrix<T>::update_matrix(const InvGammaMatrix<T>& inv_gamma, const SPLINE_G0_TYPE& spline_G0) {
  const int nop = inv_gamma.matrix().num_rows();
  const int N = matrix_.num_cols();

  sanity_check(spline_G0);

  if (nop==0) return;

  invA0.resize_values_not_retained(nop, N);
  G0_left.resize_values_not_retained(N, nop);
  G0_inv_gamma.resize_values_not_retained(N, nop);

  pl.resize(nop);
  for (int l=0; l<nop; ++l) {
    pl[l] = inv_gamma.pos_in_invA(l);
  }
  for (int l=0; l<nop; ++l) {
    for (int j=0; j<N; ++j) {
      invA0(l, j) = matrix_(pl[l],j);
    }
  }
  for (int l=0; l<nop; ++l) {
    for (int i = 0; i < N; ++i) {
      G0_left(i, l) = eval_Gij(*this, spline_G0, i, pl[l]);
    }
  }

  mygemm(-1.0, G0_left, inv_gamma.matrix(), 0.0, G0_inv_gamma);
  mygemm(1.0, G0_inv_gamma, invA0, 1.0, matrix_);

  coeff_A.resize(N);
  std::fill(coeff_A.begin(), coeff_A.end(), 1.0);

  //std::cout << "debug Ngamma " << nop << std::endl;
  for (int l=0; l < nop; ++l) {
    //std::cout << " iop " << l << " " << inv_gamma.alpha(l) << " " <<  inv_gamma.alpha0(l) << std::endl;
    const int i_row = pl[l];
    assert(inv_gamma.alpha0(l)==alpha_at(i_row));
    coeff_A[i_row] = 1.0/(
                             1.0 +
                                 gamma_func(
                                     eval_f(inv_gamma.alpha(l)),
                                     eval_f(inv_gamma.alpha0(l))
                                 )
    );
  }

  //could be vectorized
  for (int i_col=0; i_col<N; ++i_col) {
    for (int i_row=0; i_row<N; ++i_row) {
      matrix_(i_row, i_col) *= coeff_A[i_row];
    }
  }

  //update alpha
  for (int l=0; l<nop; ++l) {
    assert(pl[l]>=0 && pl[l]<N);
    alpha_[pl[l]] = inv_gamma.alpha(l);
  }
}

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
InvAMatrixFlavors<T>::remove_rows_cols(const std::vector<double>& times) {
  for (int flavor=0; flavor<sub_matrices_.size(); ++flavor) {
    sub_matrices_[flavor].remove_rows_cols(times);
  }
}

template<typename T>
template<typename SPLINE_G0_TYPE>
void
InvAMatrixFlavors<T>::sanity_check(const SPLINE_G0_TYPE& spline_G0, general_U_matrix<T>* p_Uijkl,
                                   const itime_vertex_container& itime_vertices) const {
#ifndef NDEBUG
  const int n_flavors = sub_matrices_.size();
  std::vector<std::vector<creator> > creators_scr(n_flavors);
  std::vector<std::vector<annihilator> > annihilators_scr(n_flavors);
  std::vector<std::vector<T> > alpha_scr(n_flavors);

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
    }
  }

  for (int flavor=0; flavor<n_flavors; ++flavor) {
    sub_matrices_[flavor].sanity_check(spline_G0);
    assert(sub_matrices_[flavor].annihilators().size()==annihilators_scr[flavor].size());
    assert(sub_matrices_[flavor].creators().size()==creators_scr[flavor].size());

    //check operators one by one
    for (int iop=0; iop<creators_scr[flavor].size(); ++iop) {
      const int pos = sub_matrices_[flavor].find_row_col(creators_scr[flavor][iop].t().time());
      assert(sub_matrices_[flavor].creators()[pos]==creators_scr[flavor][iop]);
      assert(sub_matrices_[flavor].annihilators()[pos]==annihilators_scr[flavor][iop]);
      assert(sub_matrices_[flavor].alpha_at(pos)==alpha_scr[flavor][iop]);
    }
  }
#endif
}

template<typename T>
template<typename SPLINE_G0_TYPE>
void InvAMatrixFlavors<T>::recompute(const SPLINE_G0_TYPE& spline_G0, bool check_error) {
  for (int flavor=0; flavor<sub_matrices_.size(); ++flavor) {
    sub_matrices_[flavor].recompute(spline_G0, check_error);
  }
}
