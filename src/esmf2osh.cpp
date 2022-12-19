#include <Omega_h_library.hpp>
#include <Omega_h_esmfWrapper.h>
#include <sstream>
#include <iostream>


int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  auto coords = new double[8];
  esmfInit();
  esmfTestMesh();
  esmfGetMeshVtxCoords(coords);
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
