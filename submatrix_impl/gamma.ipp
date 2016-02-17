#include "../submatrix.hpp"

/*
 * Implementation of InvGammaMatrix<T>
 */

template<typename T>
int InvGammaMatrix<T>::find_row_col_gamma(int pos_A) const {
  bool found = false;
  for (int i=0; i<row_col_info_.size(); ++i) {
    if (pos_A==pos_in_invA(i)) {
      found = true;
      return i;
    }
  }
  assert(found);
  throw std::runtime_error("Fatal error in InvGammaMatrix<T>::find_row_col_gamma(int pos_A): something went wrong!");
}

//Using Eq. (22)
template<typename T>
template<typename SPLINE_G0_TYPE>
T InvGammaMatrix<T>::eval_Gammaij(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0, int row_gamma, int col_gamma) const {
  if (row_gamma==col_gamma) {
    T small_gamma = gamma_func(
        eval_f(alpha(row_gamma)), eval_f(alpha0(row_gamma))
    );
    T Gij = eval_Gij_gamma(invA, spline_G0, row_gamma, col_gamma);
    return Gij - (1.0+small_gamma)/small_gamma;
  } else {
    return eval_Gij_gamma(invA, spline_G0, row_gamma, col_gamma);
  }
}

template<typename T>
template<typename SPLINE_G0_TYPE>
T InvGammaMatrix<T>::eval_Gij_gamma(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0, int row_gamma, int col_gamma) const {

  assert(pos_in_invA(row_gamma)>=0);
  assert(pos_in_invA(col_gamma)>=0);
  return eval_Gij(invA, spline_G0, pos_in_invA(row_gamma), pos_in_invA(col_gamma));
}

template<typename T>
template<typename SPLINE_G0_TYPE>
bool InvGammaMatrix<T>::sanity_check(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0) const {
  bool result = true;
#ifndef NDEBUG
  const int N = row_col_info_.size();
  assert(matrix_.num_rows()==N);
  assert(matrix_.num_cols()==N);

  if (N==0) return true;

  //check
  for (int j=0; j<N; ++j) {
    const int pj = pos_in_invA(j);
    assert(alpha0(j)==invA.alpha_at(pj));
  }

  //check gamma^-1
  alps::numeric::matrix<T> gamma(N,N,(T)0.0);
  for (int j=0; j<N; ++j) {
    for (int i=0; i<N; ++i) {
      gamma(i,j) = eval_Gij_gamma(invA, spline_G0, i, j);
    }
    const T tmp = gamma_func(eval_f(alpha(j)), eval_f(alpha0(j)));
    gamma(j,j) -= (1.0+tmp)/tmp;
  }
  alps::numeric::inverse_in_place(gamma);
  assert(alps::numeric::norm_square(gamma-matrix_)/(1.0*N*N)<1E-5);
  result = result && (alps::numeric::norm_square(gamma-matrix_)/(1.0*N*N)<1E-5);
#endif
  return result;
}

template<typename T>
void InvGammaMatrix<T>::clear() {
  resize(0);
}

template<typename T>
void InvGammaMatrix<T>::swap_rows_and_cols(int i, int j) {
  assert(i<row_col_info_.size());
  assert(j<row_col_info_.size());
  std::swap(row_col_info_[i], row_col_info_[j]);
  blas_swap_rows(matrix_,i, j);
  blas_swap_cols(matrix_,i, j);
}

template<typename T>
void InvGammaMatrix<T>::resize(size_t new_size) {
  assert(new_size>=0);
  assert(matrix_.num_rows()==matrix_.num_cols());
  assert(matrix_.num_rows()==matrix_.num_cols());

  matrix_.resize(new_size, new_size);
  row_col_info_.resize(new_size);
}

//trys to add rows and cols to Gamma
//returns the determinant ratio of A.
template<typename T>
template<typename SPLINE_G0_TYPE>
T InvGammaMatrix<T>::try_add(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0, const std::vector<OperatorToBeUpdated<T> >& ops_ins) {
  assert(matrix_.num_rows()==matrix_.num_cols());
  const int nop = matrix_.num_cols();
  const int nop_add = ops_ins.size();

  T gamma_prod = 1.0;
  for (int iop=0; iop<nop_add; ++iop) {
    //assert(invA.find_row_col(ops_ins[iop].op_t_)==ops_ins[iop].pos_in_A_);
    assert(ops_ins[iop].alpha0_==ops_ins[iop].alpha_current_);
    gamma_prod *= -gamma_func<T>(eval_f(ops_ins[iop].alpha_new_), eval_f(ops_ins[iop].alpha0_));
    row_col_info_.push_back(boost::make_tuple(ops_ins[iop].pos_in_A_, ops_ins[iop].alpha0_, ops_ins[iop].alpha_new_));
  }

  std::vector<int> rows_in_A(nop), rows_in_A2(nop_add);
  for(unsigned int i=0;i<nop;++i) {
    rows_in_A[i] = pos_in_invA(i);
  }
  for(unsigned int i=0;i<nop_add;++i) {
    rows_in_A2[i] = pos_in_invA(i+nop);
  }

  G_n_n.resize_values_not_retained(nop_add, nop_add);
  G_n_j.resize_values_not_retained(nop_add, nop);
  G_j_n.resize_values_not_retained(nop, nop_add);
  for(unsigned int i=0;i<nop;++i) {
    alps::numeric::submatrix_view<T> view(G_n_j, 0, i, nop_add, 1);
    invA.eval_Gij_col_part(spline_G0, rows_in_A2, pos_in_invA(i), view);
  }
  if (nop>0) {
    for (size_t iv=0; iv<nop_add; ++iv) {
      alps::numeric::submatrix_view<T> view(G_j_n, 0, iv, nop, 1);
      invA.eval_Gij_col_part(spline_G0, rows_in_A, pos_in_invA(iv+nop), view);
    }
  }
  for (size_t iv2=0; iv2<nop_add; ++iv2) {
    alps::numeric::submatrix_view<T> view(G_n_n, 0, iv2, nop_add, 1);
    invA.eval_Gij_col_part(spline_G0, rows_in_A2, pos_in_invA(iv2+nop), view);
    T small_gamma = gamma_func(eval_f(alpha(iv2+nop)), eval_f(alpha0(iv2+nop)));
    G_n_n(iv2, iv2) -= (1.0+small_gamma)/small_gamma;
  }

  return gamma_prod*compute_det_ratio_up(G_j_n, G_n_j, G_n_n, matrix_);
}

template<typename T>
void InvGammaMatrix<T>::perform_add() {
  //note: matrix_ is resized and then updated.
  compute_inverse_matrix_up2(G_j_n, G_n_j, G_n_n, matrix_, matrix_);
}

template<typename T>
void InvGammaMatrix<T>::reject_add() {
  row_col_info_.resize(matrix_.num_cols());
}

//it is quit unlikely that one inserts an operator which was removed already before.
//So, we assume alpha0==ALPHA_NON_INT, alpha_current!=ALPHA_NON_INT, alpha_new==ALPHA_NON_INT;
template<typename T>
template<typename SPLINE_G0_TYPE>
T InvGammaMatrix<T>::try_remove(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_g0, const std::vector<OperatorToBeUpdated<T> >& ops_rem) {
  const int nop_rem = ops_rem.size();

  T gamma_prod = 1.0;
  rows_cols_removed.resize(nop_rem);
  for (int iop=0; iop<nop_rem; ++iop) {
    assert(ops_rem[iop].alpha0_ == ALPHA_NON_INT);
    assert(ops_rem[iop].alpha_current_ != ALPHA_NON_INT);
    assert(ops_rem[iop].alpha_new_ == ALPHA_NON_INT);

    rows_cols_removed[iop] = find_row_col_gamma(ops_rem[iop].pos_in_A_);
    gamma_prod *= -gamma_func<T>(eval_f(ops_rem[iop].alpha_current_), eval_f(ops_rem[iop].alpha0_));
  }

  std::sort(rows_cols_removed.begin(), rows_cols_removed.end());
  return compute_det_ratio_down(nop_rem, rows_cols_removed, matrix_)/gamma_prod;
}

template<typename T>
void InvGammaMatrix<T>::perform_remove() {
  const int nop = matrix_.num_rows();
  const int nop_rem = rows_cols_removed.size();

  //update gamma^{-1}
  std::vector<std::pair<int,int> > rows_cols_swap_list;
  compute_inverse_matrix_down(nop_rem, rows_cols_removed, matrix_, rows_cols_swap_list);

  //remove operators
  for (int swap=0; swap<rows_cols_swap_list.size(); ++swap) {
    std::swap(row_col_info_[rows_cols_swap_list[swap].first], row_col_info_[rows_cols_swap_list[swap].second]);
  }
  row_col_info_.resize(nop-nop_rem);
}

template<typename T>
void InvGammaMatrix<T>::reject_remove() {
  //nothing to do
}

template<typename T>
template<typename SPLINE_G0_TYPE>
T InvGammaMatrix<T>::try_add_remove(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0,
                                    const std::vector<OperatorToBeUpdated<T> >& ops_ins,
                                    const std::vector<OperatorToBeUpdated<T> >& ops_rem) {
  const int nop = matrix_.num_rows();
  const int nop_add = ops_ins.size();
  const int nop_rem = ops_rem.size();
  const int nop_unchanged = nop-nop_rem;

  //move rows and cols to be removed to the end.
  //this will change the positions of the remaining rows and cols.
  T gamma_prod_rem = 1.0;
  rows_cols_removed.resize(nop_rem);
  for (int iop=0; iop<nop_rem; ++iop) {
    assert(ops_rem[iop].alpha_new_ == ALPHA_NON_INT);
    rows_cols_removed[iop] = find_row_col_gamma(ops_rem[iop].pos_in_A_);
    gamma_prod_rem *= -gamma_func<T>(eval_f(ops_rem[iop].alpha_current_), eval_f(ops_rem[iop].alpha0_));
  }
  std::sort(rows_cols_removed.begin(), rows_cols_removed.end());
  for (int swap=0; swap<nop_rem; ++swap) {
    swap_rows_and_cols(rows_cols_removed[nop_rem-1-swap], nop-1-swap);
  }

  T gamma_prod_add = 1.0;
  for (int iop=0; iop<nop_add; ++iop) {
    gamma_prod_add *= -gamma_func<T>(eval_f(ops_ins[iop].alpha_new_), eval_f(ops_ins[iop].alpha0_));
    row_col_info_.push_back(boost::make_tuple(ops_ins[iop].pos_in_A_, ops_ins[iop].alpha0_, ops_ins[iop].alpha_new_));
  }

  //At this point, row_col_info_ looks like the following:
  // op_stay_unchanged1, op_stay_unchanged2, ..., op_to_be_removed1, op_to_be_removed2, ..., op_to_be_added1, ...

  //compute the values of new entries
  G_n_n.resize_values_not_retained(nop_add, nop_add);
  G_n_j.resize_values_not_retained(nop_add, nop_unchanged);
  G_j_n.resize_values_not_retained(nop_unchanged, nop_add);
  for(unsigned int i=0;i<nop_unchanged;++i) {
    for (size_t iv=0; iv<nop_add; ++iv) {
      G_n_j(iv,i) = mycast<T>(eval_Gammaij(invA, spline_G0, nop+iv, i));
    }
  }
  for(unsigned int i=0;i<nop_unchanged;++i){
    for (size_t iv=0; iv<nop_add; ++iv) {
      G_j_n(i,iv) = mycast<T>(eval_Gammaij(invA, spline_G0, i, nop+iv));
    }
  }
  for (size_t iv2=0; iv2<nop_add; ++iv2) {
    for (size_t iv=0; iv<nop_add; ++iv) {
      G_n_n(iv, iv2) = mycast<T>(eval_Gammaij(invA, spline_G0, nop+iv, nop+iv2));
    }
  }

  return (gamma_prod_add/gamma_prod_rem)*
      compute_det_ratio_replace_rows_cols(matrix_, G_j_n, G_n_j, G_n_n, Mmat, inv_tSp);
}

template<typename T>
void InvGammaMatrix<T>::perform_add_remove() {
  const int nop = matrix_.num_rows();
  const int nop_add = G_n_n.num_cols();
  const int nop_unchanged = G_j_n.num_rows();
  const int nop_new = nop_unchanged+nop_add;

  compute_inverse_matrix_replace_rows_cols(matrix_, G_j_n, G_n_j, G_n_n, Mmat, inv_tSp);
  assert(matrix_.num_cols()==nop_new);
  std::vector<row_col_info_type> row_col_info_new(nop_new);
  for (int iop=0; iop<nop_unchanged; ++iop) {
    row_col_info_new[iop] = row_col_info_[iop];
  }
  for (int iop=0; iop<nop_add; ++iop) {
    row_col_info_new[iop+nop_unchanged] = row_col_info_[nop+iop];
  }
  std::swap(row_col_info_new, row_col_info_);
}

template<typename T>
void InvGammaMatrix<T>::reject_add_remove() {
  row_col_info_.resize(matrix_.num_cols());
}
