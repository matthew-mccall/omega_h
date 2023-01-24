#include "Omega_h_library.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_array.hpp"
#include "Omega_h_array_ops.hpp"

using namespace Omega_h;

static void test_fan_and_funnel() {
  OMEGA_H_CHECK(invert_funnel(LOs({0, 0, 1, 1, 2, 2}), 3) == LOs({0, 2, 4, 6}));
  OMEGA_H_CHECK(invert_fan(LOs({0, 2, 4, 6})) == LOs({0, 0, 1, 1, 2, 2}));
  OMEGA_H_CHECK(invert_funnel(LOs({0, 0, 0, 2, 2, 2}), 3) == LOs({0, 3, 3, 6}));
  OMEGA_H_CHECK(invert_fan(LOs({0, 3, 3, 6})) == LOs({0, 0, 0, 2, 2, 2}));
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(std::string(lib.version()) == OMEGA_H_SEMVER);
  test_fan_and_funnel();
  return 0;
}
