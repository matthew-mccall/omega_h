#include "Omega_h_library.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_array.hpp"
#include "Omega_h_array_ops.hpp"

using namespace Omega_h;

static void testFan() {
  const auto a = invert_fan(LOs({0, 2, 4, 6}));
  OMEGA_H_CHECK(a == LOs({0, 0, 1, 1, 2, 2}));
  const auto b = invert_fan(LOs({0, 2, 4, 6}));
  OMEGA_H_CHECK(b == LOs({0, 0, 1, 1, 2, 2}));
  const auto c = invert_fan(LOs({0, 2, 4, 6}));
  OMEGA_H_CHECK(c == LOs({0, 0, 1, 1, 2, 2}));
  const auto d = invert_fan(LOs({0, 2, 4, 6})); //sometimes hangs here
  OMEGA_H_CHECK(d == LOs({0, 0, 1, 1, 2, 2}));
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(std::string(lib.version()) == OMEGA_H_SEMVER);
  testFan();
  return 0;
}
