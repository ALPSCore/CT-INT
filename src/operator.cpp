//
// Created by H. Shinaoka on 2015/12/23.
//
#include "operator.hpp"

bool operator==(const creator& op1, const creator& op2) {
  return (op1.s()==op2.s()) && (op1.t()==op2.t()) && (op1.flavor()==op2.flavor());
}

bool operator==(const annihilator& op1, const annihilator& op2) {
  return (op1.s()==op2.s()) && (op1.t()==op2.t()) && (op1.flavor()==op2.flavor());
}
