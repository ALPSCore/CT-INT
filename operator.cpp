//
// Created by H. Shinaoka on 2015/12/23.
//
#include "operator.hpp"

frequency_t c_or_cdagger::nm_;
bool c_or_cdagger::use_static_exp_;
unsigned int c_or_cdagger::ntau_;
double c_or_cdagger::beta_;
double *c_or_cdagger::omegan_;
std::complex<double> *c_or_cdagger::exp_iomegan_tau_;

void c_or_cdagger::initialize_simulation(const alps::params &p)
{
  beta_=p["BETA"];
  nm_=p["NMATSUBARA_MEASUREMENTS"] | (p["NMATSUBARA"]|p["N_MATSUBARA"]);
  omegan_ = new double[nm_];
  for(unsigned int i=0;i<nm_;++i) {
    omegan_[i]=(2.*i+1.)*M_PI/beta_;
  }
  if(p.defined("TAU_DISCRETIZATION_FOR_EXP")) {
    ntau_=p["TAU_DISCRETIZATION_FOR_EXP"];
    use_static_exp_=true;
    exp_iomegan_tau_=new std::complex<double> [2*nm_*ntau_];
    if(exp_iomegan_tau_==0){throw std::runtime_error("not enough memory for computing exp!"); }
    std::cout<<"starting computation of exp values for measurement"<<std::endl;
    for(unsigned int i=0;i<ntau_;++i){
      double tau=i*beta_/(double)ntau_;
      for(unsigned int o=0;o<nm_;++o)
        exp_iomegan_tau_[2*nm_*i + o] = std::complex<double>(cos(omegan_[o]*tau), sin(omegan_[o]*tau));
      for(unsigned int o=0;o<nm_;++o)
        exp_iomegan_tau_[2*nm_*i + nm_ + o] = std::complex<double>(cos(omegan_[o]*tau), -sin(omegan_[o]*tau));
    }
    std::cout<<"done exp computation."<<std::endl;
  } else {
    use_static_exp_=false;
  }
}

bool operator==(const creator& op1, const creator& op2) {
  return (op1.s()==op2.s()) && (op1.t()==op2.t()) && (op1.flavor()==op2.flavor());
}

bool operator==(const annihilator& op1, const annihilator& op2) {
  return (op1.s()==op2.s()) && (op1.t()==op2.t()) && (op1.flavor()==op2.flavor());
}
