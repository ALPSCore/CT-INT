#pragma once

#include <algorithm>
#include <fstream>
#include "boost/multi_array.hpp"
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/cstdint.hpp>

#include "types.h"
#include "util.h"
#include <alps/params.hpp>


namespace alps {
    namespace ctint {

        typedef size_t vertex_t;
        typedef size_t af_t;

        class itime_vertex;
        class all_type;
        class density_type;
        class non_density_type;
        class non_density_type_in_window;

        template<class T>
        class vertex_definition
        {
        public:
            vertex_definition() : rank_(-1), num_af_states_(-1), Uval_(0.0), id_(-1) {
            }

            vertex_definition(size_t rank, size_t num_af_states, std::vector<spin_t>& flavors, std::vector<size_t>& sites, T Uval, boost::multi_array<T,2>& alpha_af_rank, int id)
              : rank_(rank),
                flavors_(flavors),
                sites_(sites),
                num_af_states_(num_af_states),
                Uval_(Uval),
                alpha_af_rank_(alpha_af_rank),
                id_(id) {
              assert(flavors_.size()==rank);
              assert(sites.size()==2*rank);
              assert(alpha_af_rank_.shape()[0]==num_af_states_);
              assert(alpha_af_rank_.shape()[1]==rank);
            };

            const std::vector<spin_t>& flavors() const {
              return flavors_;
            };

            const std::vector<size_t>& sites() const {
              return sites_;
            };

            T Uval() const {
              return Uval_;
            }

            size_t rank() const {
              return rank_;
            }

            size_t num_af_states() const {
              return num_af_states_;
            }

            T get_alpha(size_t af_state, size_t idx_rank) const {
              return alpha_af_rank_[af_state][idx_rank];
            }

            /**
             * This returns true if site(2*i)=site(2*i+1) for all i
             */
            bool is_density_type() const {
              bool flag = true;
              for (int i_rank=0; i_rank<rank_; ++i_rank) {
                if (sites_[2*i_rank] != sites_[2*i_rank+1]) {
                  flag = false;
                  break;
                }
              }
              return flag;
            }

            bool is_truely_non_density_type() const {
              bool flag = true;
              for (int i_rank=0; i_rank<rank_; ++i_rank) {
                flag = flag && (sites_[2*i_rank] != sites_[2*i_rank+1]);
              }
              return flag;
            }

            int id() const {return id_;}

        private:
            size_t rank_;
            std::vector<spin_t> flavors_;
            std::vector<size_t> sites_;
            size_t num_af_states_;
            T Uval_;
            boost::multi_array<T,2> alpha_af_rank_;//first index addresses af spin state, second one addresses (cdagger c)
            int id_;
        };


        template<class T>
        std::ostream& operator<<(std::ostream& os, const vertex_definition<T>& v) {
          os << " rank= " << v.rank();
          os << " flavors = ";
          for (auto f : v.flavors()) {
            os << f;
          }
          os << " sites = ";
          for (auto f : v.sites()) {
            os << f;
          }
          os << " num_af_states = " << v.num_af_states();
          for (int iaf = 0; iaf < v.num_af_states(); ++iaf) {
            for (int rank = 0; rank < v.rank(); ++rank) {
              std::cout << " ( " << iaf << "," << rank << "," << v.get_alpha(iaf, rank) << ") ";
            }
          }

          return os;
        }


//Data structure for general two-body interactions for a multi-orbital cluster impurity problem
        template<class T>
        class general_U_matrix {
        public:
            typedef T value_type;

            general_U_matrix(const alps::params &parms) :
              ns_(parms["model.sites"]),
              nf_(parms["model.spins"]) {
              if (parms.supplied("model.U_matrix_file") && !parms.supplied("model.U")) {
                load_vertex_from_file(parms);
              } else if (!parms.supplied("model.U_matrix_file") && parms.supplied("model.U")) {
                set_onsite_U(parms["model.U"]);
              } else {
                throw std::runtime_error("Do not set model.U_matrix_file and model.U at the same time!");
              }
            }

            void load_vertex_from_file(const alps::params &parms) {
              if (!parms.supplied("model.U_matrix_file")) {
                throw std::runtime_error("Error: model.U_matrix_file is not set!");
              }
              std::string ufilename(parms["model.U_matrix_file"].template as<std::string>());
              std::ifstream ifs(ufilename.c_str());
              if (!ifs.is_open()) {
                throw std::runtime_error(ufilename+" does not exist!");
              }
              ifs >> num_nonzero_;

              //temporary
              size_t rank, num_af_states;
              T Uval_;
              std::complex<double> Uval_cmplx;
              std::vector<size_t> site_indices_;//site indices (ijkl)
              std::vector<spin_t> flavor_indices_;//flavor indices for c^dagger c
              boost::multi_array<T,2> alpha_;//the first index is auxially spin, the second denotes (ij) or (kl).
              boost::multi_array<std::complex<double>,2> alpha_cmplx;//tempolary

              for (unsigned int idx=0; idx<num_nonzero_; ++idx) {
                size_t itmp;
                ifs >> itmp >> rank >> num_af_states >> Uval_cmplx;
                Uval_ = mycast<T>(Uval_cmplx);
                assert(itmp==idx);
                assert(rank==2);

                site_indices_.resize(2*rank);
                flavor_indices_.resize(rank);
                alpha_.resize(boost::extents[num_af_states][rank]);
                alpha_cmplx.resize(boost::extents[num_af_states][rank]);

                for (size_t i_op=0; i_op<2*rank; ++i_op) {
                  ifs >> site_indices_[i_op];//i, j, k, l
                  if (site_indices_[i_op] < 0 || site_indices_[i_op] >= ns_) {
                    throw std::runtime_error("Wrong site index is given in the definition of vertex.");
                  }
                }
                for (size_t i_rank=0; i_rank<rank; ++i_rank) {
                  ifs >> flavor_indices_[i_rank];
                  if (flavor_indices_[i_rank] < 0 || flavor_indices_[i_rank] >= nf_) {
                    throw std::runtime_error("Wrong flavor index is given in the definition of vertex.");
                  }
                }
                for (size_t i_rank=0; i_rank<rank; ++i_rank) {
                  for (size_t iaf=0; iaf<num_af_states; ++iaf) {
                    ifs >> alpha_cmplx[iaf][i_rank];
                  }
                }
                for (size_t i_rank=0; i_rank<rank; ++i_rank) {
                  for (size_t iaf = 0; iaf < num_af_states; ++iaf) {
                    alpha_[iaf][i_rank] = mycast<T>(static_cast<std::complex<double> >(alpha_cmplx[iaf][i_rank]));
                  }
                }

                vertex_list.push_back(vertex_definition<T>(rank, num_af_states, flavor_indices_, site_indices_, Uval_, alpha_, idx));
              }

              find_non_density_vertices();
            }

            //for unit test. Single-band Hubbard model (supported only the cases with flavor = 2)
            general_U_matrix(int n_site, double U, double delta = 1e-2) : ns_(n_site), nf_(2), num_nonzero_(n_site) {
              set_onsite_U(U, delta);
            }

            void set_onsite_U(double U, double delta = 1e-2) {
              num_nonzero_ = ns_;

              if (nf_ != 1 && nf_ != 2) {
                throw std::runtime_error("model.spins must be either 1 or 2 when model.U is specified.");
              }

              const int num_af_states = 2;
              const int rank = 2;

              std::vector<size_t> site_indices_;//site indices (ijkl)
              std::vector<spin_t> flavor_indices_;//flavor indices for c^dagger c
              boost::multi_array<T,2> alpha_;//the first index is auxially spin, the second denotes (ij) or (kl).
              boost::multi_array<std::complex<double>,2> alpha_cmplx;//tempolary

              site_indices_.resize(2*rank);
              flavor_indices_.resize(rank);
              alpha_.resize(boost::extents[num_af_states][rank]);

              auto num_physical_sites = nf_ == 2 ? ns_ : ns_/2;
              for (int site=0; site<num_physical_sites; ++site) {
                if (nf_ == 1) {
                  for (size_t i_op=0; i_op<rank; ++i_op) {
                    site_indices_[i_op] = 2*site;
                  }
                  for (size_t i_op=rank; i_op<2*rank; ++i_op) {
                    site_indices_[i_op] = 2*site + 1;
                  }
                  flavor_indices_[0] = 0;
                  flavor_indices_[1] = 0;

                  //up af spin
                  alpha_[0][0] = 1+delta;
                  alpha_[0][1] = -delta;

                  //down af spin
                  alpha_[1][0] = -delta;
                  alpha_[1][1] = 1+delta;
                } else {
                  for (size_t i_op=0; i_op<2*rank; ++i_op) {
                    site_indices_[i_op] = site;
                  }
                  flavor_indices_[0] = 0;
                  flavor_indices_[1] = 1;

                  //up af spin
                  alpha_[0][0] = 1+delta;
                  alpha_[0][1] = -delta;

                  //down af spin
                  alpha_[1][0] = -delta;
                  alpha_[1][1] = 1+delta;
                }

                vertex_list.push_back(vertex_definition<T>(rank, num_af_states, flavor_indices_, site_indices_, U, alpha_, site));
              }

              find_non_density_vertices();

            }

            size_t n_vertex_type() const{return vertex_list.size();}
            spin_t nf()const {return nf_;}
            spin_t ns()const {return ns_;}

            vertex_definition<T>& get_vertex(size_t vertex_idx) {
              assert(vertex_idx<n_vertex_type());
              return vertex_list[vertex_idx];
            }

            const vertex_definition<T>& get_vertex(size_t vertex_idx) const {
              assert(vertex_idx<n_vertex_type());
              return vertex_list[vertex_idx];
            }

            const std::vector<vertex_definition<T> >& get_vertices() const {
              return vertex_list;
            }

            const std::vector<vertex_definition<T> >& get_vertices(all_type& pred) const {
              return vertex_list;
            }

            const std::vector<vertex_definition<T> >& get_vertices(non_density_type& pred) const {
              return non_density_vertices;
            }

            const std::vector<vertex_definition<T> >& get_vertices(density_type& pred) const {
              return density_vertices;
            }

            const std::vector<vertex_definition<T> >& get_vertices(non_density_type_in_window& pred) const {
              return non_density_vertices;
            }

            const std::vector<vertex_definition<T> >& get_non_density_vertex_defs() const {
              return non_density_vertices;
            }

            int num_vertex_type(all_type& pred) const {
              return n_vertex_type();
            }

            int num_vertex_type(non_density_type& pred) const {
              return non_density_vertices.size();
            }

            int num_vertex_type(density_type& pred) const {
              return density_vertices.size();
            }

            int num_vertex_type(non_density_type_in_window& pred) const {
              return non_density_vertices.size();
            }

            int num_density_vertex_type() const {
              return density_vertices.size();
            }

            const std::vector<bool>& get_is_truely_non_density_type() const {
              return is_truely_non_density_type;
            }

        private:
            unsigned int ns_, nf_, num_nonzero_;
            std::vector<vertex_definition<T> > vertex_list, non_density_vertices, density_vertices;
            //std::vector<int> non_density_vertices;
            std::vector<bool> is_density_type, is_truely_non_density_type;

            void find_non_density_vertices() {
              is_density_type.resize(vertex_list.size());
              is_truely_non_density_type.resize(vertex_list.size());
              non_density_vertices.clear();
              for (int iv=0; iv<vertex_list.size(); ++iv) {
                is_density_type[iv] = vertex_list[iv].is_density_type();
                is_truely_non_density_type[iv] = vertex_list[iv].is_truely_non_density_type();
                if (!is_density_type[iv]) {
                  non_density_vertices.push_back(vertex_list[iv]);
                } else {
                  density_vertices.push_back(vertex_list[iv]);
                }
              }
            }
        };

//to remember what vertices are on the imaginary time axis..
        class itime_vertex {
        public:
            itime_vertex()
              : vertex_type_(-1),
                af_state_(-1),
                rank_(-1),
                time_(-1),
                is_density_type_(false),
                is_non_interacting_(false),
                unique_id_(0)
            {}

            itime_vertex(int vertex_type, int af_state, double time, int rank, bool is_density_type)
              : vertex_type_(vertex_type),
                af_state_(af_state),
                rank_(rank),
                time_(time),
                is_density_type_(is_density_type),
                is_non_interacting_(false),
                unique_id_(0)
            {}

            int af_state() const { return af_state_; }
            int vertex_type() const {return vertex_type_;}
            int type() const {return vertex_type_;}
            int rank() const {return rank_;}
            double time() const {return time_;}
            void set_time(double new_time) {time_ = new_time;}
            void set_af_state(int new_af_state) {af_state_ = new_af_state;}
            bool is_density_type() const {return is_density_type_;}
            bool is_non_interacting() const { return is_non_interacting_;}
            void set_non_interacting() { is_non_interacting_ = true;}
            void set_interacting() { is_non_interacting_ = false;}
            void set_unique_id(my_uint64 id) {unique_id_ = id;}
            void set_vertex_type(int vertex_type) {vertex_type_ = vertex_type;}
            my_uint64 unique_id() const {return unique_id_;}

        private:
            int vertex_type_, af_state_, rank_;
            double time_;
            bool is_density_type_; //, is_truely_non_density_type_;
            bool is_non_interacting_;
            my_uint64 unique_id_;
        };

        template<class V>
        class ItimeVertexContainer : private std::vector<V> {
        public:
            ItimeVertexContainer() : std::vector<V>() {};
            ItimeVertexContainer(const ItimeVertexContainer& x) : std::vector<V>(x) {};
            ItimeVertexContainer(size_t count) : std::vector<V>(count) {};

            ItimeVertexContainer<V>& operator=(const ItimeVertexContainer<V>& x) {
              std::vector<V>::operator=(x);
              return *this;
            }

            typedef typename std::vector<V>::iterator iterator;
            typedef typename std::vector<V>::const_iterator const_iterator;

            int num_interacting() const {
              int count = 0;
              for (int iv=0; iv<size(); ++iv) {
                if(!operator[](iv).is_non_interacting()) ++count;
              }
              return count;
            }

            //this is to be killed.
            using std::vector<V>::size;

            using std::vector<V>::operator[];
            //using std::vector<V>::push_back;
            using std::vector<V>::resize;
            using std::vector<V>::pop_back;
            using std::vector<V>::begin;
            using std::vector<V>::end;
            using std::vector<V>::erase;

            //const V& operator[](unsigned int index) const {
            //return operator[](index);
            //}

            void push_back(const V& x) {
              bool equal_time = false;
              for (int iv=0; iv<size(); ++iv) {
                if (operator[](iv).time()==x.time()) {
                  equal_time = true;
                  break;
                }
              }
              if (equal_time) {
                std::runtime_error("ItimeVertexContainer::push_back: you try to insert a vertex at the imaginary time where there is already a vertex.");
              }
              std::vector<V>::push_back(x);
            }
        };

        typedef ItimeVertexContainer<itime_vertex> itime_vertex_container;

        inline bool operator<(const itime_vertex& v1, const itime_vertex& v2) {
          return (v1.time()<v2.time());
        }
//template<class T>
//itime_vertex generate_itime_vertex(vertex_definition<T> vertex_df
//template<class V>
//int num_non_denisty_vertices(const V& itime_vertices) {
        //int r = 0;
        //const std::vector<bool>& tmp = Uijkl.get_is_denisty_type();
        //for (typename V::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
        //if (tmp[it->type()]) ++r;
        //}
        //return r;
//}

        class density_type : public std::unary_function<itime_vertex,bool> {
        public:
            bool operator()(itime_vertex v) {return v.is_density_type();}

            double random_time(double random01, double beta) const {
              assert(random01>=0 && random01<=1);
              return random01*beta;
            }
        };

        class non_density_type : public std::unary_function<itime_vertex,bool> {
        public:
            bool operator()(itime_vertex v) {return !v.is_density_type();}

            double random_time(double random01, double beta) const {
              assert(random01>=0 && random01<=1);
              return random01*beta;
            }
        };

//this does not include correlated hopping.
/*
class truely_non_density_type : public std::unary_function<itime_vertex,bool> {
public:
    bool operator()(itime_vertex v) {
      if (v.is_density_type())
        return false;

      bool flag = true;
      for (int ir=0; ir<v.rank(); ++ir) {
        flag = flag && v.
      }
      return v.is_truely_non_density_type();
    }

    double random_time(double random01, double beta) const {
      assert(random01>=0 && random01<=1);
      return random01*beta;
    }
};
*/


        class non_density_type_in_window : public std::unary_function<itime_vertex,bool> {
        public:
            non_density_type_in_window() :
              ts_(0),
              w_(0),
              t_small1_(0),
              t_large1_(1E+100),
              t_small2_(0),
              t_large2_(1E+100) {}

            non_density_type_in_window(double ts, double w, double beta) : ts_(ts), w_(w), beta_(beta) {
              assert(w<=beta);
              if (ts+w<=beta) {
                t_small1_ = t_small2_= ts;
                t_large1_ = t_large2_ = ts+w;
              } else {
                t_small1_ = ts;
                t_large1_ = beta;
                t_small2_ = 0.0;
                t_large2_ = w+ts-beta;
              }
              //std::cout << "t_s1,t_l1" << t_small1_ << " " << t_large1_ << std::endl;
              //std::cout << "t_s2,t_l2" << t_small2_ << " " << t_large2_ << std::endl;
            }

            bool operator()(itime_vertex v) const {
              const double t = v.time();
              return !v.is_density_type() && ((t_small1_<=t && t<=t_large1_) || (t_small2_<=t && t<=t_large2_));
            }

            double random_time(double random01, double beta)  const{
              assert(random01>=0 && random01<=1);
              assert(beta_==beta);
              assert(w_>0);
              double t = random01*w_+ts_;
              if (t>beta_) t -= beta_;
              return t;
            }

            double width() {
              assert(w_>0);
              return w_;
            }

        private:
            double ts_, w_, beta_;
            double t_small1_, t_large1_;
            double t_small2_, t_large2_;
        };

        class all_type : public std::unary_function<itime_vertex,bool> {
        public:
            bool operator()(itime_vertex v) {return true;}

            double random_time(double random01, double beta)  const{
              assert(random01>=0 && random01<=1);
              return random01*beta;
            }
        };

        template<class T, class R, class UnaryPredicate>
        std::vector<itime_vertex> generate_valid_vertex_pair(const general_U_matrix<T>& Uijkl, const std::vector<std::tuple<int,int> >& valid_pair,
                                                             R& random01, double beta, const UnaryPredicate& pred) {
          const int n_vertices_add = 2;

          std::vector<itime_vertex> itime_vertices;
          itime_vertices.reserve(n_vertices_add);

          if (valid_pair.size()==0)
            return std::vector<itime_vertex>();

          int v1,v2;
          std::tie(v1,v2) = valid_pair[static_cast<int>(random01()*valid_pair.size())];
          int vtypes[] = {v1, v2};

          for (int iv=0; iv<n_vertices_add; ++iv) {
            const double time = pred.random_time(random01(), beta);
            const int v_type = vtypes[iv];
            const int rank = Uijkl.get_vertices()[v_type].rank();
            const int af_state = random01()*Uijkl.get_vertices()[v_type].num_af_states();
            itime_vertices.push_back(itime_vertex(v_type, af_state, time, rank, false));
            if(Uijkl.get_vertices()[v_type].is_density_type()) {
              throw std::logic_error("Error found density type vertex");
            }
          }
          return itime_vertices;
        }

        template<class T, class R, class P>
        std::vector<itime_vertex> generate_valid_vertex_pair2(const general_U_matrix<T>& Uijkl, const std::pair<int,int> v_pair,
                                                              R& random01, double beta, const P& normalized_prob_dist) {
          const int n_vertices_add = 2;

          std::vector<itime_vertex> itime_vertices;
          itime_vertices.reserve(n_vertices_add);

          int vtypes[] = {v_pair.first, v_pair.second};
          double times[2];
          times[0] = mymod(beta*random01(), beta);
          times[1] = mymod(times[0]+gen_rand_rejection_method(normalized_prob_dist, normalized_prob_dist(0.0), random01, beta), beta);
          assert(vtypes[0]<vtypes[1]);

          for (int iv=0; iv<n_vertices_add; ++iv) {
            const double time = times[iv];
            const int v_type = vtypes[iv];
            const int rank = Uijkl.get_vertices()[v_type].rank();
            const int af_state = random01()*Uijkl.get_vertices()[v_type].num_af_states();
            itime_vertices.push_back(itime_vertex(v_type, af_state, time, rank, false));
            if(Uijkl.get_vertices()[v_type].is_density_type()) {
              throw std::logic_error("Error found density type vertex");
            }
          }
          return itime_vertices;
        }


        template<class T, class R, class UnaryPredicate>
        std::vector<itime_vertex> generate_itime_vertices(const general_U_matrix<T>& Uijkl, R& random01, double beta, int n_vertices_add, UnaryPredicate pred) {
          std::vector<itime_vertex> itime_vertices;
          itime_vertices.reserve(n_vertices_add);

          const std::vector<vertex_definition<T> >& valid_vs = Uijkl.get_vertices(pred);
          const int n_valid_vs = valid_vs.size();
          if (n_valid_vs==0)
            return std::vector<itime_vertex>();

          for (int iv=0; iv<n_vertices_add; ++iv) {
            const double time = pred.random_time(random01(), beta);
            const int iv_rnd = static_cast<int>(random01()*n_valid_vs);
            const int v_type = valid_vs[iv_rnd].id();
            const int rank = valid_vs[iv_rnd].rank();
            const int af_state = static_cast<size_t>(random01()*valid_vs[iv_rnd].num_af_states());
            itime_vertices.push_back(itime_vertex(v_type, af_state, time, rank, valid_vs[iv_rnd].is_density_type()));
          }
          return itime_vertices;
        }

        template<class UnaryPredicate>
        int count_valid_vertex_pair(const std::vector<itime_vertex>& itime_vertices, const boost::multi_array<bool,2>& valid_pair_flag, const UnaryPredicate& window) {
          std::vector<int> vs_win; vs_win.reserve(itime_vertices.size());
          std::vector<int> pos; pos.reserve(itime_vertices.size());

          int idx = 0;
          for (std::vector<itime_vertex>::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
            if (window(*it) && !it->is_density_type()) {
              vs_win.push_back(it->type());
              pos[vs_win.size()-1] = idx;
            }
            ++idx;
          }

          const int N = vs_win.size();
          int N_valid_pair = 0;
          for (int iv=0; iv<N; ++iv) {
            for (int iv2=0; iv2<N; ++iv2) {
              if (iv<=iv2)
                continue;

              if (valid_pair_flag[vs_win[iv]][vs_win[iv2]])
                ++N_valid_pair;
            }
          }

          return N_valid_pair;
        }

        template<class UnaryPredicate, class R>
        std::vector<int>
        pick_up_valid_vertex_pair(const std::vector<itime_vertex>& itime_vertices, const boost::multi_array<bool,2>& valid_pair_flag, const UnaryPredicate& window, R& random01, int& N_valid_pair) {
          std::vector<int> vs_win; vs_win.reserve(itime_vertices.size());
          std::vector<int> pos; pos.reserve(itime_vertices.size());

          int idx = 0;
          for (std::vector<itime_vertex>::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
            if (window(*it) && !it->is_density_type()) {
              vs_win.push_back(it->type());
              pos.push_back(idx);
            }
            ++idx;
          }

          const int N = vs_win.size();
          std::vector<std::pair<int,int> > pos_pair;
          N_valid_pair = 0;
          for (int iv=0; iv<N; ++iv) {
            for (int iv2=0; iv2<N; ++iv2) {
              if (iv<=iv2)
                continue;

              if (valid_pair_flag[vs_win[iv]][vs_win[iv2]]) {
                ++N_valid_pair;
                pos_pair.push_back(std::make_pair(pos[iv],pos[iv2]));
              }
            }
          }

          if (N_valid_pair==0)
            return std::vector<int>();

          assert(pos_pair.size()==N_valid_pair);
          //const int selected = static_cast<int>(random01()*N_valid_pair);
          //std::cout << " debug " << selected << " " << N_valid_pair << " " << pos_pair.size() << std::endl;
          //assert(selected<N_valid_pair);
          std::pair<int,int> tmp = pos_pair[static_cast<int>(random01()*N_valid_pair)];
          assert(!itime_vertices[tmp.first].is_density_type());
          assert(!itime_vertices[tmp.second].is_density_type());
          std::vector<int> r(2);
          r[0] = tmp.first;
          r[1] = tmp.second;
          return r;
        }

        template<class P, class R>
        std::pair<int,int>
        pick_up_valid_vertex_pair2(const itime_vertex_container& itime_vertices, std::pair<int,int> v_pair,
                                   double beta, P& p, R& random01, int& N_valid_pair, double& F) {

          std::vector<itime_vertex> v1, v2;
          std::vector<int> pos_v1, pos_v2;
          std::vector<double> prob;

          const int Nv = itime_vertices.size();
          for (int iv=0; iv<Nv; ++iv) {
            if (itime_vertices[iv].is_non_interacting()) {
              continue;
            }
            if (itime_vertices[iv].type()==v_pair.first) {
              v1.push_back(itime_vertices[iv]);
              pos_v1.push_back(iv);
            } else if (itime_vertices[iv].type()==v_pair.second) {
              v2.push_back(itime_vertices[iv]);
              pos_v2.push_back(iv);
            }
          }

          N_valid_pair = v1.size()*v2.size();
          if (N_valid_pair==0) {
            return std::make_pair(-1,-1);
          }

          assert(v1.size()==pos_v1.size());
          assert(v2.size()==pos_v2.size());
          prob.resize(0); prob.reserve(v1.size()*v2.size());
          for (int iv1=0; iv1<v1.size(); ++iv1) {
            for (int iv2=0; iv2<v2.size(); ++iv2) {
              assert(itime_vertices[pos_v2[iv2]].type()==v_pair.second);
              assert(itime_vertices[pos_v1[iv1]].type()==v_pair.first);
              prob.push_back(p.bare_value(
                mymod(itime_vertices[pos_v2[iv2]].time()-itime_vertices[pos_v1[iv1]].time(), beta)
              ));
            }
          }

          F = std::accumulate(prob.begin(), prob.end(), 0.0);
          const int idx = boost::random::discrete_distribution<>(prob.begin(), prob.end())(random01.engine());
          assert(idx/v2.size()<v1.size());
          assert(idx%v2.size()<v2.size());
          return std::make_pair(pos_v1[idx/v2.size()], pos_v2[idx%v2.size()]);
        }


        template<class R, class UnaryPredicate>
        std::vector<int> pick_up_itime_vertices(const itime_vertex_container& itime_vertices,
                                                R& random01,
                                                int n_vertices_rem,
                                                UnaryPredicate pred) {
          const int n_active_vertices = std::count_if(itime_vertices.begin(), itime_vertices.end(), pred);
          if (n_active_vertices<n_vertices_rem)
            return std::vector<int>();

          std::vector<int> pos(n_active_vertices);
          int idx=0;
          for (int iv=0; iv<itime_vertices.size(); ++iv) {
            if (pred(itime_vertices[iv])) {
              assert(idx<pos.size());
              pos[idx] = iv;
              ++idx;
            }
          }
          assert(idx==n_active_vertices);

          const std::vector<int>& indices = pickup_a_few_numbers(n_active_vertices, n_vertices_rem, random01);
          std::vector<int> indices2(n_vertices_rem);
          for (int i=0; i<n_vertices_rem; ++i) {
            assert(indices[i]<pos.size());
            indices2[i] = pos[indices[i]];
          }

#ifndef NDEBUG
          for (int i=0; i<n_vertices_rem; ++i) {
            assert(indices2[i]<itime_vertices.size());
            assert(pred(itime_vertices[indices2[i]]));
          }
#endif

          return indices2;
        };

//void
//find_valid_pair_multi_vertex_update(const std::vector<std::vector<std::valarray<int> > >& quantum_numbers, std::vector<boost::tuple<int,int,int,int> >& v_pair, boost::multi_array<bool,4>& v_pair_flag);
        template<class T>
        void
        find_valid_pair_multi_vertex_update(const std::vector<vertex_definition<T> >& vertex_defs, const std::vector<std::vector<std::valarray<int> > >& quantum_numbers, std::vector<std::pair<int,int> >& v_pair, boost::multi_array<bool,2>& v_pair_flag) {
          v_pair.resize(0);
          v_pair_flag.resize(boost::extents[quantum_numbers.size()][quantum_numbers.size()]);

          std::valarray<int> qn(quantum_numbers[0][0].size());
          const int i_af1=0, i_af2=0;
          for (int v1=0; v1<quantum_numbers.size(); ++v1) {
            if (is_all_zero<int>(quantum_numbers[v1][i_af1]))
              continue;
            for (int v2=0; v2<quantum_numbers.size(); ++v2) {
              if (is_all_zero<int>(quantum_numbers[v2][i_af2]))
                continue;
              qn = quantum_numbers[v1][i_af1] + quantum_numbers[v2][i_af2];
              bool flag = is_all_zero<int>(qn);

              v_pair_flag[v1][v2] = flag;
              if (v1 < v2 && flag)
                v_pair.push_back(std::make_pair(v1, v2));
            }
          }
        };

        template<class V>
        void print_vertices(std::ostream &os, const V &v) {
          os << std::endl;
          for (int iv = 0; iv < v.size(); ++iv) {
            os << " iv = " << iv;
            os << " type= " << v[iv].type();
            os << " rank= " << v[iv].rank();
            os << " af_state= " << v[iv].af_state();
            os << " time= " << v[iv].time();
            os << std::endl;
          }
        }

        template<typename T>
        void load_config(std::ifstream &is, const general_U_matrix<T>& Uijkl, itime_vertex_container &itime_vertices) {
          itime_vertices.resize(0);

          int Nv;
          is >> Nv;
          if (Nv==0) return;

          for (int iv=0; iv<Nv; ++iv) {
            int iv_in, type, af_state;
            double time;
            is >> iv_in >> type >> af_state >> time;
            //std::cout << "type, af_state, time " << type << " " << af_state << " " << time << std::endl;
            const vertex_definition<T>& vdef = Uijkl.get_vertex(type);
            itime_vertices.push_back(itime_vertex(type, af_state, time, vdef.rank(), vdef.is_density_type()));
          }
        }

        inline
        std::ostream &operator<<(std::ostream &os, const itime_vertex &v) {
          os << " type= " << v.type();
          os << " rank= " << v.rank();
          os << " af_state= " << v.af_state();
          os << " time= " << v.time();
          return os;
        }

        inline
        void dump(std::ostream &os, const itime_vertex_container &itime_vertices) {
          const int Nv = itime_vertices.size();
          std::cout << "Nv = " << Nv << std::endl;
          os << Nv << std::endl;
          int iv = 0;
          for (itime_vertex_container::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
            os << iv << " " << it->type() << " " << it->af_state() << " " << it->time() << std::endl;
            ++iv;
          }
        }

    }
}
