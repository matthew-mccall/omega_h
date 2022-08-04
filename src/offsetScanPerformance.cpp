#include "Omega_h_array.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_profile.hpp"

using namespace Omega_h;

static void test_scan(size_t size) {
  auto w = Write<LO>(size, 1, "foo");
  auto r = read(w);
  std::stringstream ss;
  ss << "offsetScan_" << size;
  std::string name = ss.str();
  ScopedTimer scanTime(name.c_str());
  offset_scan(r, name);
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(std::string(lib.version()) == OMEGA_H_SEMVER);
  auto size = std::size_t(1) << 21;
  ScopedTimer totalTime("total time");
  for(int i=0; i<16; i++) {
    size = size >> 1;
    printf("%d size %lu\n", i, size);
    for(int j=0; j<100; j++) {
      test_scan(size);
    }
  }
}
