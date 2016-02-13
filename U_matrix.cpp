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

void print_vertices(std::ostream &os, const std::vector<itime_vertex> &v) {
  os << std::endl;
  for (int iv=0; iv<v.size(); ++iv) {
    os << " iv = " << iv;
    os << " type= " << v[iv].type();
    os << " rank= " << v[iv].rank();
    os << " af_state= " << v[iv].af_state();
    os << " time= " << v[iv].time();
    os << std::endl;
  }
}


void dump(std::ostream &os, const itime_vertex_container &itime_vertices) {
  const int Nv = itime_vertices.size();
  std::cout << "Nv = " << Nv << std::endl;
  int iv = 0;
  for (itime_vertex_container::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
    os << iv << " " << it->type() << " " << it->af_state() << " " << it->time() << std::endl;
    ++iv;
  }
}


