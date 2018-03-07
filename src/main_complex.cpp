#include "main.hpp"

int main(int argc, char** argv) {
  using namespace alps::ctint;
  return run_simulation<InteractionExpansion<complex_number_solver> >(argc, argv);
}

