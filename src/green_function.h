/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Matthias Troyer <troyer@comp-phys.org>
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H
#include "types.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include <algorithm>

//#include <alps/mcbase.hpp>
#include "util.h"

#include <boost/format.hpp>


#ifdef USE_MPI
#include <mpi.h>
#endif 

//#include <alps/hdf5/archive.hpp>
//#include <alps/hdf5/pointer.hpp>
//#include <alps/hdf5/complex.hpp>

#include "U_matrix.h"

//forward declaration
//class general_U_matrix;
//class vertex_definition;

//types
typedef std::valarray<int> quantum_number_t;

//Matsubara GF: use T=std::complex<double>
//Imaginary time: use T=double
template <typename T> class green_function{

public:
  //construction and destruction, assignement and copy constructor
  ///constructor: how many time slices, how many sites, how many flavors
  green_function(unsigned int ntime, unsigned int nsite, unsigned int nflavor):nt_(ntime), ns_(nsite), nf_(nflavor),
  ntnsns_(ntime*nsite*nsite), ntns_(ntime*nsite){
    val_=new T[nt_*ns_*ns_*nf_];
    err_=new T[nt_*ns_*ns_*nf_];
  }
  ///specialization: constructor for problems with only one site
  green_function(unsigned int ntime, unsigned int nflavor):nt_(ntime), ns_(1), nf_(nflavor),
  ntnsns_(ntime), ntns_(ntime){
    val_=new T[nt_*nf_];
    err_=new T[nt_*nf_];
  }
  ///destructor
  ~green_function(){
    delete [] val_;
    delete [] err_;
  }
  ///copy constructor
  green_function(const green_function &g):nt_(g.nt_), ns_(g.ns_), nf_(g.nf_), ntnsns_(g.ntnsns_), ntns_(g.ntns_){
    val_=new T[nt_*ns_*ns_*nf_];
    err_=new T[nt_*ns_*ns_*nf_];
    operator=(g);
  }
  ///operator= (assignement operator)
  const green_function &operator=(const green_function &g){
    memcpy(val_, g(), sizeof(T)*nt_*ns_*ns_*nf_);
    memcpy(err_, g.error(), sizeof(T)*nt_*ns_*ns_*nf_);
    return *this;
  }
  void clear(){ memset(val_, 0, ns_*ns_*nt_*nf_*sizeof(T)); }
  //access of vectors and elements
  ///specialization for only one site: access element with given time and flavor
  inline T &operator()(unsigned int t, unsigned int flavor){return val_[t+nt_*flavor];}
  ///specialization for only one site: return const reference to element with given time and flavor
  inline const T &operator()(unsigned int t, unsigned int flavor)const{return val_[t+nt_*flavor];}

  ///return an entire vector of times for a given flavor
  inline T *operator()(unsigned int flavor){return val_+ntnsns_*flavor;}

  //error access
  inline T &error(unsigned int t, unsigned int flavor){return err_[t+nt_*flavor];}
  inline const T &error(unsigned int t, unsigned int flavor)const{return err_[t+nt_*flavor];}
  inline T *errors(unsigned int flavor){return err_+nt_*flavor;}
  ///access element with given time, site 1, site 2, and flavor
  inline T &operator()(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor){return val_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
  ///access element with given time, site 1, site 2, and flavor (const reference)
  inline const T &operator()(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor)const{return val_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
  ///return an entire vector of imaginary time values for a given site 1, site2, flavor
  inline T *operator()(unsigned int site1, unsigned int site2, unsigned int flavor){return val_+nt_*site1+ntns_*site2+ntnsns_*flavor;}

  inline T &error(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor){return err_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
  inline const T &error(unsigned int t, unsigned int site1, unsigned int site2, unsigned int flavor)const{return err_[t+nt_*site1+ntns_*site2+ntnsns_*flavor];}
  inline T *errors(unsigned int site1, unsigned int site2, unsigned int flavor){return err_+nt_*site1+ntns_*site2+ntnsns_*flavor;}

  ///get all values at once
  inline const T *operator()() const {return val_;}
  ///get all errors at once
  inline const T *error() const {return err_;}

  //size information
  ///how many flavors do we have? (flavors are usually spins, GF of different flavors are zero)
  unsigned int nflavor()const{return nf_;}
  ///return # of sites
  unsigned int nsite()const{return ns_;}
  ///return # of imaginary time values
  unsigned int ntime()const{return nt_;}
  ///return # of matsubara frequencies. Exactly equivalent to ntime().
  ///In the case of a Matsubara GF 'ntime' sounds odd -> define 'nfreq' instead.
  unsigned int nfreq()const{return nt_;} //nfreq is an alias to ntime - more intuitive use for Matsubara GF
  void read(const char *filename);
  void write(const char *filename) const;
  void write_hdf5(alps::hdf5::archive &ar, const std::string &path) const{
    ar<<alps::make_pvp(path+"/nt",nt_);
    ar<<alps::make_pvp(path+"/ns",ns_);
    ar<<alps::make_pvp(path+"/nf",nf_);
    std::stringstream subpath; subpath<<path<<"/values/mean";
    ar<<alps::make_pvp(subpath.str(), val_, nt_*ns_*ns_*nf_); // the nondiagonal components are needed for realspace representation of multisite problems
  }
  void read_hdf5(alps::hdf5::archive &ar, const std::string &path) {
    unsigned int nt, ns, nf;
    clear();
    ar>>alps::make_pvp(path+"/nt",nt);
    ar>>alps::make_pvp(path+"/ns",ns);
    ar>>alps::make_pvp(path+"/nf",nf);
    if(nt!=nt_ || ns!=ns_ || nf!=nf_){ std::cerr<<path<<" nt: "<<nt_<<" new: "<<nt<<" ns: "<<ns_<<" "<<ns<<" nf: "<<nf_<<" "<<nf<<" dimensions do not match."<<std::endl; throw std::runtime_error("Green's function read in: dimensions do not match."); }
    /*
    if (ns==1) {
      for(unsigned int i=0;i<nf_;++i){
        std::stringstream subpath; subpath<<path<<"/"<<i<<"/mean/value";
        ar>>alps::make_pvp(subpath.str(), val_+i*nt_, nt_);
        //currently we're not writing the error.
        //std::stringstream subpath_e; subpath_e<<path<<"/"<<i<<"/mean/error";
        //ar<<alps::make_pvp(subpath_e.str(), err_+i*nt_, nt_);
      }
    } else {*/
    std::stringstream subpath; subpath<<path<<"/values/mean";
    ar>>alps::make_pvp(subpath.str(), val_, nt_*ns_*ns_*nf_);
    //}
//    std::cerr << "3";
  }

  //returns if G(site1, site2) is zero.
  bool is_zero(size_t site1, size_t site2, spin_t flavor, double eps) const {
    assert(site1<nsite() && site2<nsite());
    bool flag = true;
    for (size_t i=0; i<nfreq(); ++i) {
      if (std::abs(operator()(i,site1,site2,flavor))>eps) {
        flag = false;
        break;
      }
    }
    return flag;
  }

  std::pair<std::vector<T>,std::vector<T> > to_multiple_vector() const;
  void from_multiple_vector(const std::pair<std::vector<T>,std::vector<T> > &mv);
#ifdef USE_MPI
  void broadcast(){
    MPI_Bcast( val_, ntnsns_*nf_*sizeof(T)/sizeof(double), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast( err_, ntnsns_*nf_*sizeof(T)/sizeof(double), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
#endif

private:
  //const values
  const unsigned int nt_; ///imag time points
  const unsigned int ns_; ///number of sites
  const unsigned int nf_; ///number of flavors
  const unsigned int ntnsns_; ///nt*ns*ns
  const unsigned int ntns_; ///nt*ns
  // the actual values and errors.
  T *val_;
  T *err_;
};

//diagonal in flavor space
//typedef green_function<std::complex<double> > matsubara_green_function_t;
//typedef green_function<double> itime_green_function_t;

//has off-diagonal elements in flavor space
//typedef full_green_function<std::complex<double> > matsubara_full_green_function_t;
//typedef full_green_function<double> itime_full_green_function_t;

///write out imag time Green function
std::ostream &operator<<(std::ostream &os, const green_function<double> &v);
///read in imag time Green function
std::istream &operator>>(std::istream &is, green_function<double> &v);
///write out Matsubara Green function
std::ostream &operator<<(std::ostream &os, const green_function<std::complex<double> > &v);
///read in Matsubara Green function
std::istream &operator>>(std::istream &is, green_function<std::complex<double> > &v);

///compute kinetic energy
double kinetic_energy(const multiple_vector_type &G_tau, const double &beta, const double &t);

template<typename T> void green_function<T>::read(const char *filename){
  std::ifstream in_file(filename);
  assert(in_file.is_open()); //this is not good enough when people use -DNDEBUG.
  if(!in_file.is_open()){ throw(std::invalid_argument("input file could not be opened!")); }
  double ignored=0;
  for(unsigned int i=0;i<nt_;++i){
    in_file>>ignored; //read first entry, which could be # matsubara frequencies, or tau-point, or N/beta*tau, or...
    for(unsigned int s0=0; s0<ns_; ++s0) {
      for(unsigned int s1=0; s1<ns_; ++s1){
        for(unsigned int f=0; f<nf_; ++f) {
          in_file>>operator()(i, s0, s1, f)>>std::ws; //read the actual value
        }
      }
    }
  }
}

template<typename T> void green_function<T>::write(const char *filename) const{
  std::ofstream out_file(filename);
  assert(out_file.is_open()); //this is not good enough when people use -DNDEBUG.
  if(!out_file.is_open()){ std::cerr<<"output file "<<filename<<" could not be opened!"<<std::endl; exit(1);}
  for(unsigned int i=0;i<nt_;++i){
    out_file << i << " ";
    for(unsigned int s0=0; s0<ns_; ++s0) {
      for(unsigned int s1=0; s1<ns_; ++s1){
        for(unsigned int f=0; f<nf_; ++f) {
          out_file<<operator()(i, s0, s1, f) << " "; 
        }
      }
    }
    out_file << std::endl;
  }
}

///for the transition period from multiple vectors to this data structure only.
template<typename T> std::pair<std::vector<T>,std::vector<T> > green_function<T>::to_multiple_vector() const{
  assert(ns_==1 && nf_<=2);
  std::pair<std::vector<T>,std::vector<T> > mv;
  mv.first.resize(nt_);
  mv.second.resize(nt_);
  for(unsigned int i=0;i<nt_;++i){
    mv.first[i]=operator()(i,0,0,0);
    mv.second[i]=nf_==1?operator()(i,0,0,0):operator()(i,0,0,1);
  }
  return mv;
}
template<typename T> void green_function<T>::from_multiple_vector(const std::pair<std::vector<T>,std::vector<T> >&mv){
  assert(ns_==1 && nf_<=2);
  for(unsigned int i=0;i<nt_;++i){
    operator()(i,0,0,0)=mv.first[i];
    if(nf_==2)
      operator()(i,0,0,1)=mv.second[i];
  }
}


enum shape_t {diagonal, blockdiagonal, nondiagonal};


/*
void print_all_green_functions(std::string const &basename, const int iteration_ctr, const matsubara_green_function_t &G0_omega,
                               const matsubara_green_function_t &G_omega, const itime_green_function_t &G0_tau, 
                               const itime_green_function_t &G_tau, const double beta, const shape_t shape=diagonal,
                               const std::string suffix="");
void print_real_green_matsubara(std::ostream &os, const matsubara_green_function_t &v, const double beta, const shape_t shape=diagonal);
void print_imag_green_matsubara(std::ostream &os, const matsubara_green_function_t &v, const double beta, const shape_t shape=diagonal);
void print_tau_green_functions(std::string const &basename, const int iteration_ctr, const itime_green_function_t &G0_tau, const itime_green_function_t &G_tau, const double beta,
                               const shape_t shape=nondiagonal, const std::string suffix="");
void print_dressed_tau_green_functions(std::string const &basename, const int iteration_ctr, const itime_green_function_t &G_tau, const double beta, 
                                       const shape_t shape=nondiagonal, const std::string suffix="");
                                       */


template<typename T>
std::pair<green_function<std::complex<double> >, green_function<T> >
read_bare_green_functions(const alps::params &parms) {
  const unsigned int n_matsubara = parms["N_MATSUBARA"];
  const unsigned int n_tau=parms["N_TAU"];
  const spin_t n_flavors(parms["FLAVORS"]);
  const unsigned int n_site(parms["SITES"]);

  green_function<std::complex<double> > bare_green_matsubara(n_matsubara, n_site, n_flavors);
  green_function<T> bare_green_itime(n_tau+1, n_site, n_flavors);

  //Load G0 in Matsubara frequency
  {
    if (parms["G0_OMEGA"] == "") {
      throw std::runtime_error("Please set G0_OMEGA");
    }
    std::string myfile(parms["G0_OMEGA"].template as<std::string>());
    std::ifstream ifs(myfile.c_str());
    if (!ifs.is_open()) {
      throw std::runtime_error(myfile+" does not exist!");
    }
    int flavor_tmp, itmp, itmp2, itmp3;
    double re, imag;
    int line = 0;
    for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
      for (std::size_t site1=0; site1<n_site; ++site1) {
        for (std::size_t site2=0; site2<n_site; ++site2) {
          for (std::size_t iomega = 0; iomega < n_matsubara; iomega++) {
            ifs >> flavor_tmp >> itmp >> itmp2 >> itmp3 >> re >> imag;
            if (flavor_tmp!=flavor) {
              throw std::runtime_error((boost::format("Ilegal format of G0_OMEGA: We expect %1% at the first column of the line %2%") % flavor % line).str().c_str());
            }
            if (itmp!=site1) {
              throw std::runtime_error((boost::format("Ilegal format of G0_OMEGA: We expect %1% at the second column of the line %2%") % site1 % line).str().c_str());
            }
            if (itmp2!=site2) {
              throw std::runtime_error((boost::format("Ilegal format of G0_OMEGA: We expect %1% at the third column of the line %2%") % site2 % line).str().c_str());
            }
            if (itmp3!=iomega) {
              throw std::runtime_error((boost::format("Ilegal format of G0_OMEGA: We expect %1% at the fourth column of the line %2%") % iomega % line).str().c_str());
            }
            bare_green_matsubara(iomega, site1, site2, flavor) = std::complex<double>(re, imag);
            ++line;
          }
        }
      }
    }
  }

  //Load G0 in imaginary time
  {
    if (parms["G0_TAU"] == "") {
      throw std::runtime_error("Please set G0_TAU");
    }
    std::string myfile(parms["G0_TAU"].template as<std::string>());
    std::ifstream ifs(myfile.c_str());
    if (!ifs.is_open()) {
      throw std::runtime_error(myfile+" does not exist!");
    }
    int flavor_tmp, itmp, itmp2, itmp3;
    double re, im;
    int line = 0;
    for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
      for (std::size_t site1=0; site1<n_site; ++site1) {
        for (std::size_t site2=0; site2<n_site; ++site2) {
          for (std::size_t itau = 0; itau < n_tau+1; itau++) {
            ifs >> flavor_tmp >> itmp >> itmp2 >> itmp3 >> re >> im;

            if (flavor_tmp != flavor) {
              throw std::runtime_error(
                  (boost::format("Ilegal format of G0_TAU: We expect %1% at the first column of the line %2%") % flavor %
                   line).str().c_str());
            }
            if (itmp != site1) {
              throw std::runtime_error(
                  (boost::format("Ilegal format of G0_TAU: We expect %1% at the second column of the line %2%") % site1 %
                   line).str().c_str());
            }
            if (itmp2 != site2) {
              throw std::runtime_error(
                  (boost::format("Ilegal format of G0_TAU: We expect %1% at the third column of the line %2%") %
                   site2 % line).str().c_str());
            }
            if (itmp3 != itau) {
              throw std::runtime_error(
                  (boost::format("Ilegal format of G0_TAU: We expect %1% at the fourth column of the line %2%") % itau %
                   line).str().c_str());
            }
            bare_green_itime(itau, site1, site2, flavor) = mycast<T>(std::complex<double>(re, im));
            ++line;
          }
        }
      }
    }
  }

  return std::pair< green_function<std::complex<double> >, green_function<T> >(bare_green_matsubara, bare_green_itime);
};


//groups(groups, sites belonging to groups)
template<class T>
void
make_groups(size_t N, const T& connected, std::vector<std::vector<size_t> >& groups, std::vector<int>& map) {
  map.resize(N);
  std::fill(map.begin(),map.end(),-1);
  groups.clear();

  for (size_t site1=0; site1<N; ++site1) {
    int connected_to_where = -1;
    for (size_t site2=0; site2<site1; ++site2) {
      if (connected[site1][site2] && map[site2]!=-1) {
        connected_to_where = map[site2];
        break;
      }
    }
    if (connected_to_where==-1) {
      //create a new group
      groups.resize(groups.size()+1);
      groups[groups.size()-1].push_back(site1);
      map[site1] = groups.size()-1;
    } else {
      groups[connected_to_where].push_back(site1);
      map[site1] = connected_to_where;
    }
  }
#ifndef NDEBUG
  {
    size_t t_sum = 0;
    for (size_t ig=0; ig<groups.size(); ++ig) {
      t_sum += groups[ig].size();
    }
    assert(t_sum==N);
    for (size_t i=0; i<N; ++i) {
      assert(map[i]>=0);
    }
  }
#endif
}

//template<class T>
//bool compare_array_size(std::vector<T>& array1, std::vector<T>& array2) {
  //return array1.size()<array2.size();
//}

template<class T>
void print_group(const std::vector<std::vector<T> >& group) {
  for (int i=0; i<group.size(); ++i) {
    std::cout << "   Group "<<i<<" consists of site ";
    for (int j=0; j<group[i].size(); ++j) {
      std::cout << group[i][j];
      if (j<group[i].size()-1)
        std::cout << ", ";
    }
    std::cout << "." << std::endl;
  }
}

template<typename T, typename S, typename G>//Expected T,S=double or T=std::complex<double>
std::vector<std::vector<quantum_number_t> >
make_quantum_numbers(const G& gf, const std::vector<vertex_definition<S> >& vertices,
                     std::vector<std::vector<std::vector<size_t> > >& groups,
                     std::vector<std::vector<int> >& group_map,
                     double eps=1E-10) {
  const size_t n_site = gf.nsite();
  const size_t n_flavors = gf.nflavor();

  //See if two sites are connected by nonzero G
  boost::multi_array<bool,2> connected(boost::extents[n_site][n_site]);
  groups.clear(); groups.resize(n_flavors);
  group_map.clear(); group_map.resize(n_flavors);
  std::vector<size_t> num_groups(n_flavors);

  //figure out how sites are connected by GF
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    for (size_t site1 = 0; site1 < n_site; ++site1) {
      for (size_t site2 = 0; site2 < n_site; ++site2) {
        connected[site1][site2] = !gf.is_zero(site1, site2, flavor, eps);
      }
    }
    make_groups(n_site, connected, groups[flavor], group_map[flavor]);
    num_groups[flavor] = groups[flavor].size();
  }

  //determine the dimension of quantum numbers
  const size_t n_dim = *std::max_element(num_groups.begin(), num_groups.end());

  //compute quantum number for each vertex
  const size_t Nv = vertices.size();
  std::vector<std::vector<quantum_number_t> > qn_vertices(Nv);
  for (size_t iv=0; iv<Nv; ++iv) {
    const vertex_definition<S>& vd = vertices[iv];
    const int num_af = vd.num_af_states();
    for (int i_af=0; i_af<num_af; ++i_af) {
      std::valarray<int> qn_diff(0, n_dim*n_flavors);
      assert(qn_diff.size()==n_dim*n_flavors);
      for (size_t i_rank=0; i_rank<vd.rank(); ++i_rank) {
        const spin_t flavor = vd.flavors()[i_rank];
        const size_t site1 = vd.sites()[2*i_rank];//c_dagger
        const size_t site2 = vd.sites()[2*i_rank+1];//c

        int PH;
        if (site1==site2) {
          //density type
          PH = std::abs(vd.get_alpha(i_af,i_rank))<std::abs(vd.get_alpha(i_af,i_rank)-1.0) ? 1 : -1;
        } else {
          //non-density type
          PH = 1;
        }

        //C^dagger
        assert(n_dim*flavor+group_map[flavor][site1]<qn_diff.size());
        qn_diff[n_dim*flavor+group_map[flavor][site1]] += PH;
        //c
        assert(n_dim*flavor+group_map[flavor][site2]<qn_diff.size());
        qn_diff[n_dim*flavor+group_map[flavor][site2]] -= PH;
      }
      qn_vertices[iv].push_back(qn_diff);
    }
  }

  return qn_vertices;
}

#endif
