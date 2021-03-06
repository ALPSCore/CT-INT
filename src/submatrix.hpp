//
// Created by H. Shinaoka on 2016/01/26.
//

#pragma once

#include <algorithm>

#include <boost/tuple/tuple.hpp>

#include <boost/lambda/lambda.hpp>
#include <boost/function.hpp>

#include "operator.hpp"
#include "U_matrix.h"
#include "fastupdate_formula.h"

namespace alps {
    namespace ctint {


//The following should be moved to U_matrix.h
        const double ALPHA_NON_INT = 1E+100;
        const int NON_INT_SPIN_STATE = -1;

/*
 * Forward definition
 */
        template<class T> class InvGammaMatrix;
        template<class T> class InvAMatrix;

        template<typename T>
        inline T eval_f(T alpha) {
          return alpha==ALPHA_NON_INT ? 1.0 : alpha/(alpha-1.0);
        }

        template<typename T>
        inline T gamma_func(T f1, T f2) {
          return (f1-f2)/f2;
        }

        template<typename T, typename SPLINE_G0_TYPE>
        T eval_Gij(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0, int row_A, int col_A);

        template<class T>
        class InvAMatrix
        {
        public:
            typedef T value_type;
            typedef boost::tuple<vertex_t,size_t,my_uint64> vertex_info_type;

            InvAMatrix();

            alps::numeric::matrix<T> &matrix() { return matrix_;}
            alps::numeric::matrix<T> const &matrix() const { return matrix_;}
            alps::numeric::matrix<T> const &tildeG() const { return tildeG_;}
            std::vector<creator> &creators(){ return creators_;}
            const std::vector<creator> &creators() const{ return creators_;}
            std::vector<annihilator> &annihilators(){ return annihilators_;}
            const std::vector<annihilator> &annihilators()const{ return annihilators_;}
            T alpha_at(int pos) const {
              assert(pos>=0 && pos<creators_.size());
              return creators_[pos].s()==annihilators_[pos].s() ? alpha_[pos] : alpha_[pos];
            }
            T bare_alpha_at(int pos) const {
              assert(pos>=0 && pos<creators_.size());
              return alpha_[pos];
            }
            void set_alpha(int pos, T new_val) {
              assert(pos>=0 && pos<creators_.size());
              alpha_[pos] = new_val;
            }
            void alpha_push_back(T new_elem) {alpha_.push_back(new_elem);}
            std::vector<vertex_info_type> &vertex_info(){ return vertex_info_;}
            const std::vector<vertex_info_type>& vertex_info() const{ return vertex_info_;}

            int find_row_col(my_uint64 v_uid, int i_rank) const {
              for(std::size_t i=0; i<creators_.size(); ++i) {
                //if (time==creators_[i].t().time() && vertex_info_[i].first==type && vertex_info_[i].second==i_rank) {
                if (vertex_uid(i)==v_uid && boost::get<1>(vertex_info_[i])==i_rank) {
                  return i;
                }
              }
              return -1;
            }

            std::vector<int> find_row_col(my_uint64 v_uid) const {
              std::vector<int> pos;
              for(std::size_t i=0; i<creators_.size(); ++i) {
                if (vertex_uid(i)==v_uid) {
                  pos.push_back(i);
                }
              }
              if (pos.size()==0) {
                throw std::logic_error("No operator found in InvAMatrix::find_row_col().");
              }
              return pos;
            }

            my_uint64 vertex_uid(int pos) const {
              assert(pos>=0 && pos<vertex_info_.size());
              return boost::get<2>(vertex_info_[pos]);
            }
            template<typename SPLINE_G0_TYPE>
            bool sanity_check(const SPLINE_G0_TYPE& spline_G0) const;
            void swap_ops(size_t i1, size_t i2) {
              assert(i1>=0 && i1<creators_.size());
              assert(i2>=0 && i2<creators_.size());
              std::swap(creators_[i1], creators_[i2]);
              std::swap(annihilators_[i1], annihilators_[i2]);
              std::swap(alpha_[i1], alpha_[i2]);
              std::swap(vertex_info_[i1], vertex_info_[i2]);
            }
            void swap_rows_cols(size_t i1, size_t i2) {
              swap_ops(i1, i2);
              blas_swap_cols(matrix_, i1, i2);
              blas_swap_rows(matrix_, i1, i2);
            }
            template<class InputIterator>
            void swap_ops2(InputIterator first, InputIterator end) {
              for (InputIterator it=first; it!=end; ++it) {
                swap_ops(it->first, it->second);
              }
            }
            template<typename SPLINE_G0_TYPE>
            void update_matrix(const InvGammaMatrix<T>& inv_gamma, const SPLINE_G0_TYPE& spline_G0);
            void remove_rows_cols(const std::vector<int>& rows_cols);
            void pop_back_op() {
              creators_.pop_back();
              annihilators_.pop_back();
              alpha_.pop_back();
              vertex_info_.pop_back();
            }
            void push_back_op(const creator& cdag_op, const annihilator& c_op, T alpha, const vertex_info_type& vertex_info);
            template<typename SPLINE_G0_TYPE>
            void extend(const SPLINE_G0_TYPE& spline_G0);
            void remove_rows_cols(const std::vector<my_uint64>& v_uid);

            template<typename SPLINE_G0_TYPE>
            void update_tildeG(const SPLINE_G0_TYPE& spline_G0) {
              tildeG_.destructive_resize(matrix_.size1(), matrix_.size1());
              std::vector<int> all_rows_cols;
              for (int i=0; i < matrix_.size1(); ++i) {
                all_rows_cols.push_back(i);
              }
              eval_Gij_cols_rows(spline_G0, all_rows_cols, all_rows_cols, tildeG_);
            }

            /*recompute A^{-1} and return det(A) and det(1-F)*/
            template<typename SPLINE_G0_TYPE>
            std::pair<T,T> recompute_matrix(const SPLINE_G0_TYPE& spline_G0, bool check_error);

            T compute_f_prod() const;

            template<typename SPLINE_G0_TYPE>
            void compute_M(alps::numeric::matrix<T>& M, const SPLINE_G0_TYPE& spline_G0) const;

            template<typename M>
            void read_Gij_col(int col, M& Gij) const;

            template<typename M>
            void read_Gij_cols_rows(const std::vector<int>& rows,
                                    const std::vector<int>& cols, M& result) const;


            template<typename SPLINE_G0_TYPE, typename M>
            void eval_Gij_cols_rows(const SPLINE_G0_TYPE& spline_G0,
                               const std::vector<int>& rows,
                               const std::vector<int>& cols, M& result) const;

            template<typename SPLINE_G0_TYPE, typename M>
            void eval_Gij_col_part(const SPLINE_G0_TYPE& spline_G0, const std::vector<int>& rows, int col, M& Gij) const;

        private:
            //compute G0 (and reuse cached data)
            template<typename SPLINE_G0_TYPE>
            alps::numeric::submatrix_view<T> compute_G0_col(const SPLINE_G0_TYPE& spline_G0, int col) const;

            alps::numeric::matrix<T> matrix_;
            alps::numeric::matrix<T> tildeG_;
            std::vector<creator> creators_;         //an array of creation operators c_dagger corresponding to the row of the matrix
            std::vector<annihilator> annihilators_; //an array of to annihilation operators c corresponding to the column of the matrix
            std::vector<T> alpha_;             //an array of doubles corresponding to the alphas of Rubtsov for the c, cdaggers at the same index.
            std::vector<vertex_info_type> vertex_info_; // an array of pairs which remember from which type of vertex operators come from. (type of vertex and rank)

            //work space for update()
            alps::numeric::matrix<T> G0_left, invA0, G0_inv_gamma;
            std::vector<int> pl;

            //cache for G0 obtained by interpolation
            mutable alps::numeric::matrix<T> G0_cache;
            mutable std::vector<int> index_G0_cache;
            mutable int num_entry_G0_cache;
        };

        template<class T>
        class InvAMatrixFlavors
        {
        public:
            InvAMatrixFlavors(int n_flavors) : sub_matrices_(n_flavors) {
              assert(n_flavors>=0);
            }

            InvAMatrix<T>& operator[](size_t flavor) {
              assert(flavor<sub_matrices_.size());
              return sub_matrices_[flavor];
            }

            template<typename SPLINE_G0_TYPE>
            std::pair<T,T> init(general_U_matrix<T>* p_Uijkl, const SPLINE_G0_TYPE& spline_G0, const itime_vertex_container& itime_vertices, int begin_index);

            const InvAMatrix<T>& operator[](size_t flavor) const {
              assert(flavor<sub_matrices_.size());
              return sub_matrices_[flavor];
            }

            size_t size() const {
              return sub_matrices_.size();
            };

            template<typename SPLINE_G0_TYPE>
            bool sanity_check(const SPLINE_G0_TYPE& spline_G0, general_U_matrix<T>* p_Uijkl, const itime_vertex_container& itime_vertices) const;

            template<typename SPLINE_G0_TYPE>
            std::pair<T,T> recompute_matrix(const SPLINE_G0_TYPE& spline_G0, bool check_error);

            typename InvAMatrix<T>::value_type determinant() {
              typename InvAMatrix<T>::value_type det=1.0;
              for (spin_t flavor=0; flavor<size(); ++flavor) {
                det *= sub_matrices_[flavor].matrix().determinant();
              }
              return det;
            }

            template<typename SPLINE_G0_TYPE>
            void add_non_interacting_vertices(general_U_matrix<T>* p_Uijkl, const SPLINE_G0_TYPE& spline_G0, const itime_vertex_container& itime_vertices, int begin_index);

            template<typename SPLINE_G0_TYPE>
            void update_tildeG(const SPLINE_G0_TYPE& spline_G0) {
              for (int flavor=0; flavor<sub_matrices_.size(); ++flavor) {
                sub_matrices_[flavor].update_tildeG(spline_G0);
              }
            }

            void remove_rows_cols(const std::vector<my_uint64>& vertex_uid);

            template<typename SPLINE_G0_TYPE>
            void update_matrix(const std::vector<InvGammaMatrix<T> >& inv_gamma_flavors, const SPLINE_G0_TYPE& spline_G0);

        private:
            std::vector<InvAMatrix<T> > sub_matrices_;
        };

        template<class T>
        class OperatorToBeUpdated {
        public:
            OperatorToBeUpdated(operator_time op_t, int pos_in_A, T alpha0, T alpha_current, T alpha_new)
              : op_t_(op_t), pos_in_A_(pos_in_A), alpha0_(alpha0), alpha_current_(alpha_current), alpha_new_(alpha_new) {};
            OperatorToBeUpdated()
              : op_t_(), pos_in_A_(-1), alpha0_(-1), alpha_current_(-1), alpha_new_(-1) {};

            operator_time op_t_;
            int pos_in_A_;
            T alpha0_, alpha_current_, alpha_new_;
        };

        template<class T>
        class InvGammaMatrix {
        public:
            typedef T value_type;
            typedef boost::tuple<int,T,T> row_col_info_type;//position in A^{-1}, alpha0, alpha_current

            InvGammaMatrix() {}

            const alps::numeric::matrix<T>& matrix() const {return matrix_;}

            //adding new rows and cols
            template<typename SPLINE_G0_TYPE>
            T try_add(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_g0, const std::vector<OperatorToBeUpdated<T> >& ops_ins);
            void perform_add();
            void reject_add();

            //removing rows and cols
            template<typename SPLINE_G0_TYPE>
            T try_remove(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_g0, const std::vector<OperatorToBeUpdated<T> >& ops_rem);
            void perform_remove();
            void reject_remove();

            //removing rows and cols and then adding new rows and cols
            template<typename SPLINE_G0_TYPE>
            T try_add_remove(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_g0, const std::vector<OperatorToBeUpdated<T> >& ops_ins,
                             const std::vector<OperatorToBeUpdated<T> >& ops_rem);
            void perform_add_remove();
            void reject_add_remove();

            void clear();

            int pos_in_invA(int row_col_in_gamma) const {
              assert(row_col_in_gamma<row_col_info_.size());
              return boost::get<0>(row_col_info_[row_col_in_gamma]);
            }

            T alpha0(int row_col_in_gamma) const {
              assert(row_col_in_gamma<row_col_info_.size());
              return boost::get<1>(row_col_info_[row_col_in_gamma]);
            }

            T alpha(int row_col_in_gamma) const {
              assert(row_col_in_gamma<row_col_info_.size());
              return boost::get<2>(row_col_info_[row_col_in_gamma]);
            }

            template<typename SPLINE_G0_TYPE>
            T eval_Gammaij(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0, int row, int col) const;

            template<typename SPLINE_G0_TYPE>
            T eval_Gij_gamma(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0, int row, int col) const;

            template<typename SPLINE_G0_TYPE>
            bool sanity_check(const InvAMatrix<T>& invA, const SPLINE_G0_TYPE& spline_G0) const;


        private:
            int noperators0;
            alps::numeric::matrix<T> matrix_;//inverse of gamma matrix
            std::vector<row_col_info_type> row_col_info_;//info about rows and cols in gamma matrix
            std::vector<creator> creators_gamma_;
            std::vector<annihilator> annihilators_gamma_;

            //workspace for addition of row and col
            alps::numeric::matrix<T> G_n_n, G_n_j, G_j_n;

            //workspace for removal of row and col
            std::vector<int> rows_cols_removed;
            alps::numeric::matrix<T> Mmat, inv_tSp;

            //auxially functions for multi-vertex insertion and removal
            int find_row_col_gamma(int pos_A) const;

            void swap_rows_and_cols(int i, int j);
            void resize(size_t new_size);
        };

        template<class T>
        class SubmatrixUpdate
        {
        public:
            typedef boost::function<T(const annihilator&,const creator&)> SPLINE_G0_TYPE;

            SubmatrixUpdate(int k_ins_max, int n_flavors, SPLINE_G0_TYPE spline_G0, general_U_matrix<T>* p_Uijkl, double beta);//, const alps::params &p);

            SubmatrixUpdate(int k_ins_max, int n_flavors, SPLINE_G0_TYPE spline_G0, general_U_matrix<T>* p_Uijkl, double beta,
                            const itime_vertex_container& itime_vertices_init);//, const alps::params &p);


            InvAMatrix<T>& submatrix(size_t flavor) {
              assert(flavor<invA_.size());
              return invA_[flavor];
            }

            const InvAMatrix<T>& submatrix(size_t flavor) const {
              assert(flavor<invA_.size());
              return invA_[flavor];
            }

            const T sign() const {
              return sign_;
            }

            size_t n_flavors() const {
              return invA_.size();
            };

            int k_ins_max() const {
              return k_ins_max_;
            }

            //this should be called when there are no non-interacting vertices
            int pert_order() const {
#ifndef NDEBUG
              for (int i=0; i<itime_vertices_.size(); ++i) {
                assert(!itime_vertices_[i].is_non_interacting());
              }
#endif
              return itime_vertices_.size();
            }

            const itime_vertex_container& itime_vertices() const {
              return itime_vertices_;
            }

            /*
             * Compute M=G0^-1 from A^{-1}
             * For measuring self-energy
             */
            void compute_M(std::vector<alps::numeric::matrix<T> >& M);

            /*
             * Compute M=G0^-1 from scratch
             * Return sign of Monte Carlo weight and weight itselft.
             * This is very expensive. Only for test use.
             */
            std::pair<T,T> compute_M_from_scratch(std::vector<alps::numeric::matrix<T> >& M);

            const InvAMatrixFlavors<T>& invA() const {
              return invA_;
            }

            //The core of submatrix update
            void init_update(const std::vector<itime_vertex>& non_int_itime_vertices);

            boost::tuple<T,T,T>
            try_spin_flip(const std::vector<int>& pos, const std::vector<int>& new_spins);

            void perform_spin_flip(const std::vector<int>& pos, const std::vector<int>& new_spins);

            void reject_spin_flip();

            void finalize_update();

            /* recomputes A^{-1} to avoid numerical errors*/
            void recompute_matrix(bool check_error);

            //for debug
            bool sanity_check();

        private:
            enum SubmatrixState {READY_FOR_UPDATE=0, TRYING_SPIN_FLIP=1};
            const int k_ins_max_;
            SPLINE_G0_TYPE spline_G0_; //for interpolation of G0
            general_U_matrix<T>* p_Uijkl_;
            const double beta_;
            //const T coeff_det;

            SubmatrixState state;
            InvAMatrixFlavors<T> invA_;

            //Monte Carlo variables
            T sign_det_A_, sign_;
            std::vector<InvGammaMatrix<T> > gamma_matrices_;
            itime_vertex_container itime_vertices_;//current configuration
            itime_vertex_container itime_vertices0_;//starting configuration for submatrix update set in init_update()

            //workspace
            T det_rat_A, sign_rat;
            std::vector<std::vector<OperatorToBeUpdated<T> > > ops_rem, ops_ins, ops_replace;//operator_time and new alpha

            //parameters
            alps::params params;

            my_uint64 gen_new_vertex_id();
            my_uint64 current_vertex_id_;
        };

#include "./submatrix_impl/common.ipp"
#include "./submatrix_impl/invA.ipp"
#include "./submatrix_impl/gamma.ipp"

    }
}
