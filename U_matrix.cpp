//
// Created by H. Shinaoka on 2015/10/30.
//
#include "U_matrix.h"

std::ostream &operator<<(std::ostream &os, const itime_vertex &v) {
  os << " type= " << v.type();
  os << " rank= " << v.rank();
  os << " af_state= " << v.af_state();
  os << " time= " << v.time();
  return os;
}

