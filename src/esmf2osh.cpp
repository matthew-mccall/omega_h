#include <Omega_h_library.hpp>
#include <sstream>
#include <iostream>

extern "C" void esmfInit();
extern "C" void esmfFinalize();
extern "C" void esmfTestMeshGet(double* coords);

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  auto coords = new double[8];
  esmfInit();
  esmfTestMeshGet(coords);
  esmfFinalize();
  std::stringstream ss;
  for(int i=0; i<4; i++) {
    ss << "(" << coords[i*2] << "," << coords[i*2+1] << ") ";
  }
  ss << "\n";
  std::cout << ss.str();
  delete [] coords;
  return 0;
}
