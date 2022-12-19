#include <Omega_h_library.hpp>
#include <Omega_h_esmfWrapper.h>
#include <sstream>
#include <iostream>


int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  auto coords = new double[8];
  auto elemVerts = new int[6];
  auto vtxIds = new int[4];
  esmfInit();
  esmfTestMesh();
  esmfGetMeshVtxCoords(coords);
  esmfGetMeshVtxIds(vtxIds);
  esmfGetMeshElemVerts(elemVerts);
  esmfFinalize();
  std::stringstream ss;
  ss << "vtxIds: ";
  for(int i=0; i<4; i++) {
    ss << vtxIds[i] << ", ";
  }
  ss << "\n";
  std::cout << ss.str();
  ss.str("");
  ss << "vtxCoords: ";
  for(int i=0; i<4; i++) {
    ss << "(" << coords[i*2] << "," << coords[i*2+1] << ") ";
  }
  ss << "\n";
  std::cout << ss.str();
  ss.str("");
  ss << "elemVerts: ";
  for(int i=0; i<2; i++) {
    ss << "(" << elemVerts[i*3] << "," << elemVerts[i*3+1] << "," << elemVerts[i*3+2] << ") ";
  }
  ss << "\n";
  std::cout << ss.str();
  delete [] coords;
  return 0;
}
