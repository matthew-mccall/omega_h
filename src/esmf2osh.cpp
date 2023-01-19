#include <Omega_h_library.hpp>
#include <Omega_h_esmfWrapper.h>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <sstream>
#include <iostream>

int fileTypeStringToInt(std::string fileType) {
  assert(fileType == "scrip" || fileType == "esmf" || fileType == "ugrid");
  if(fileType == "scrip") return 2;
  if(fileType == "esmf") return 3;
  if(fileType == "ugrid") return 5;
  else return 0; //UNKNOWN
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("input.osh");
  cmdline.add_arg<std::string>("output.exo");
  auto& fileTypeFlag = cmdline.add_flag("--in-mesh-type", "scrip|esmf|ugrid");
  fileTypeFlag.add_arg<std::string>("type");
  if (!cmdline.parse(comm, &argc, argv) ||
      !Omega_h::CmdLine::check_empty(comm, argc, argv)) {
    cmdline.show_help(comm, argv);
    return -1;
  }
  auto inpath = cmdline.get<std::string>("input.osh");
  auto outpath = cmdline.get<std::string>("output.exo");
  auto fileTypeString = cmdline.get<std::string>("--in-mesh-type", "type");
  auto fileType = fileTypeStringToInt(fileTypeString);
  esmfInit();
  std::string meshFileName(argv[1]);
  esmfLoadMesh(inpath.c_str(), inpath.length(), fileType);
  int dim, numVerts, numElms;
  esmfGetMeshInfo(&dim, &numVerts, &numElms);
  std::cout << "dim, numVerts, numElms: "
            << dim << ", " << numVerts << ", " << numElms << "\n";
  esmfFinalize();
  return 0;
}
