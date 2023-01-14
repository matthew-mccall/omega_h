#include <Omega_h_library.hpp>
#include <Omega_h_esmfWrapper.h>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <sstream>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  esmfInit();
  std::string meshFileName(argv[1]);
  esmfLoadMesh(meshFileName.c_str(), meshFileName.length());
  int dim, numVerts, numElms;
  esmfGetMeshInfo(&dim, &numVerts, &numElms);
  std::cout << dim << " " << numVerts << " " << numElms << "\n";
  esmfFinalize();
  return 0;
}
