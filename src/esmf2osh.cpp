#include <Omega_h_library.hpp>
#include <Omega_h_esmfWrapper.h>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <sstream>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("input.osh");
  cmdline.add_arg<std::string>("output.exo");
  if (!cmdline.parse(comm, &argc, argv) ||
      !Omega_h::CmdLine::check_empty(comm, argc, argv)) {
    cmdline.show_help(comm, argv);
    return -1;
  }
  auto inpath = cmdline.get<std::string>("input.osh");
  auto outpath = cmdline.get<std::string>("output.exo");
  esmfInit();
  std::string meshFileName(argv[1]);
  esmfLoadMesh(inpath.c_str(), inpath.length());
  int dim, numVerts, numElms;
  esmfGetMeshInfo(&dim, &numVerts, &numElms);
  std::cout << dim << " " << numVerts << " " << numElms << "\n";
  esmfFinalize();
  return 0;
}
