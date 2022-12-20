#include <Omega_h_library.hpp>
#include <Omega_h_esmfWrapper.h>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
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

  //create element-to-vtx device array
  //the esmf 'elemVerts' array contains indices into the array of vertex ids
  //the esmf vertex ids start at one instead of zero
  Omega_h::HostWrite<Omega_h::LO> elemVerts_hw(6);
  for(int i=0; i<elemVerts_hw.size(); i++)
    elemVerts_hw[i] = (vtxIds[elemVerts[i]-1]-1);
  auto elemVerts_dr = Omega_h::read(elemVerts_hw.write());

  Omega_h::HostWrite<Omega_h::Real> vtxCoords_hw(8);
  for(int i=0; i<vtxCoords_hw.size(); i++)
    vtxCoords_hw[i] = coords[i];
  auto vtxCoords_dr = Omega_h::read(vtxCoords_hw.write());

  auto mesh = Omega_h::Mesh(&lib);
  Omega_h::build_from_elems_and_coords(&mesh, OMEGA_H_SIMPLEX, 2,
      elemVerts_dr, vtxCoords_dr);
  Omega_h::binary::write("twoTri.osh", &mesh);

  delete [] coords;
  delete [] elemVerts;
  delete [] vtxIds;
  esmfFinalize();
  return 0;
}
