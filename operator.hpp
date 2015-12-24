/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
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

#ifndef DMFT_QMC_WEAK_COUPLING_OPERATOR_H
#define DMFT_QMC_WEAK_COUPLING_OPERATOR_H

#include <alps/ngs/params.hpp>
#include "types.h"

extern "C" void vdsin_(const int *n, const double *a, double *y);
extern "C" void vdcos_(const int *n, const double *a, double *y);
extern "C" void vdsincos_(const int *n, const double *a, double *s, double *c);
extern "C" void vrda_sincos_(const int *n, const double *a, double *s, double *c);

template<class T>
class operator_time_tmpl
{
public:
    operator_time_tmpl() : time_(0), small_idx_(0) {}
    operator_time_tmpl(T time, int small_idx) : time_(time), small_idx_(small_idx) {}

    inline T time() const {return time_;}
    inline int small_index() const {return small_idx_;}

    inline void set_time(T time) { time_ = time;}
    inline void set_small_index(int small_idx) { small_idx_ = small_idx;}

private:
    T time_;
    int small_idx_;
};

typedef operator_time_tmpl<itime_t> operator_time;

template<class T>
bool operator<(const operator_time_tmpl<T>& t1, const operator_time_tmpl<T>& t2) {
  if (t1.time()==t2.time()) {
    return (t1.small_index()<t2.small_index());
  } else {
    return (t1.time()<t2.time());
  }
}

/*creation and annihilation operator class*/
class c_or_cdagger   //represents a creation operator or an annihilation operator
{ 
public:
  c_or_cdagger( const spin_t z,const site_t s, const operator_time t, const frequency_t n_matsubara)
  {
    s_ = s;
    z_ = z;
    t_ = t;
    nm_ = n_matsubara;
    exp_computed_ = false;
    exp_iomegat_ = 0;
  }
  
  
  
  ~c_or_cdagger()
  {
    if(exp_computed_)
      delete [] exp_iomegat_;
  }
  
  
  
  const c_or_cdagger & operator=(const c_or_cdagger &c)
  {
    if(this != &c){
      s_ = c.s_;
      z_ = c.z_;
      t_ = c.t_;
      if(use_static_exp_)
        exp_iomegat_=c.exp_iomegat_;
      else {
         if (exp_computed_ && c.exp_computed_) {
          memcpy(exp_iomegat_, c.exp_iomegat_, sizeof(std::complex<double>)*c.nm_);
        }
        else if (exp_computed_ && (!c.exp_computed_)) {
          delete [] exp_iomegat_;
        }
        else if ((!exp_computed_) && c.exp_computed_) {
          exp_iomegat_ = new std::complex<double>[c.nm_];
          memcpy(exp_iomegat_, c.exp_iomegat_, sizeof(std::complex<double>)*c.nm_);
        }
      }
      nm_=c.nm_;
      exp_computed_=c.exp_computed_;          
    }
    return *this;
  }
  
  
  
  c_or_cdagger(const c_or_cdagger & c)
  {
    exp_computed_=false;
    operator=(c);
  }
  


  inline const spin_t &flavor() const {return z_;}
  inline spin_t &flavor() {return z_;}
  inline const operator_time &t() const{return t_;}
  inline operator_time &t() {return t_;}
  inline const site_t &s() const {return s_;}
  inline void flavor(spin_t z){z_=z;}
  inline void s(site_t s){s_=s;}
  inline const std::complex<double> * exp_iomegat() const {return exp_iomegat_;} 
  //contains exp(iomegat) if its a creator, exp(-iomegat) if it's an annihilator.
  static void initialize_simulation(const alps::params &parms);
  
  
  static const std::complex<double> *exp_iomegan_tau(const double &tau) 
  {
    int taun=(int)(tau*ntau_/beta_); 
    return &(exp_iomegan_tau_[taun*2*nm_]);
  }
  

  static const std::complex<double> *exp_min_iomegan_tau(const double &tau) 
  {
    int taun=(int)(tau*ntau_/beta_);
    return &(exp_iomegan_tau_[taun*2*nm_ + nm_]);
  }
  
  
  
private:
  site_t s_;      //this vertex's site
  operator_time t_;     //its imaginary time point
  spin_t z_;      //its flavor or flavor
  static unsigned int nm_;        //number of matsubara frequencies
  std::complex<double> *exp_iomegat_;
  bool exp_computed_;
  static bool use_static_exp_; //do we compute the exps once only or for each time slice?
  static unsigned int ntau_;
  static double beta_;
  static double *omegan_;
  static std::complex<double> *exp_iomegan_tau_;//huge field with exp(i omega_n tau) for discretized taus.



public:
  void set_time(operator_time time) {
    t_ = time;
  }


  void compute_exp(const frequency_t n_matsubara, const int sign)
  {
    if(!use_static_exp_){
      if(!exp_computed_){
        double* sin_array = new double[n_matsubara];
        double* cos_array = new double[n_matsubara];
        //ACML vector functions
#ifdef ACML
        int one=1;
        double arg_array[n_matsubara];
        int nm=n_matsubara;
        memcpy(arg_array, omegan_, n_matsubara*sizeof(double));
        dscal_(&n_matsubara, &t_, arg_array, &one);
        vrda_sincos_(&nm, arg_array, sin_array, cos_array);
#else 
        //MKL vector functions
#ifdef MKL
        int one=1;
        double arg_array = new double[n_matsubara];
        int nm=n_matsubara;
        memcpy(arg_array, omegan_, n_matsubara*sizeof(double));
        dscal_(&n_matsubara, &t_, arg_array, &one);
        vdsincos_(&nm, arg_array, sin_array, cos_array);
        delete [] arg_array;
#else
        //NO vector functions
        for(frequency_t o=0;o<n_matsubara;++o){
          cos_array[o]=cos(omegan_[o]*t_.time());
          sin_array[o]=sin(omegan_[o]*t_.time());
        }
#endif
#endif    
        exp_iomegat_=new std::complex<double>[n_matsubara];
        for(frequency_t o=0;o<n_matsubara;++o)
          exp_iomegat_[o] = std::complex<double>(cos_array[o], sign*sin_array[o]);
        exp_computed_=true;
        delete[] sin_array;
        delete[] cos_array;
      }
    } else { //use static exp
      int taun=(int)(t_.time()*ntau_/beta_);
      assert(taun<ntau_);
      if(sign==1) 
        exp_iomegat_=&(exp_iomegan_tau_[taun*2*nm_]);
      else
        exp_iomegat_=&(exp_iomegan_tau_[taun*2*nm_ + nm_]);
    }


  }

};// creator, annihilator;

class creator : public c_or_cdagger {
public:
    using c_or_cdagger::operator=;
    creator(const spin_t z,const site_t s, const operator_time t, const frequency_t n_matsubara)
            : c_or_cdagger(z, s, t, n_matsubara){};
    ~creator(){};
};

class annihilator : public c_or_cdagger {
public:
    using c_or_cdagger::operator=;
    annihilator(const spin_t z,const site_t s, const operator_time t, const frequency_t n_matsubara)
            : c_or_cdagger(z, s, t, n_matsubara){};
    ~annihilator(){};
};

#endif

