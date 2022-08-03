#include "Omega_h_array.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_for.hpp"

using namespace Omega_h;

static void test_set() {
  auto w = Write<LO>(100, "foo");
  OMEGA_H_CHECK(w.size() == 100);
  OMEGA_H_CHECK(w.name() == "foo");
  w.set(0,42);
  OMEGA_H_CHECK(42 == w.get(0));
  w.set(88,1010);
  OMEGA_H_CHECK(1010 == w.get(88));
  auto f = OMEGA_H_LAMBDA(LO i) { 
    if( i == 88 ) assert( w[i] == 1010 );
    if( i == 0 ) assert( w[i] == 42 );
    printf("%d %d\n", i, w[i]);
  };
  parallel_for(w.size(), f, "check_vals");
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(std::string(lib.version()) == OMEGA_H_SEMVER);
  test_set();
}
