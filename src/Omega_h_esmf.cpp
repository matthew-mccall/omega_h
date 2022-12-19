#include "Omega_h_file.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cctype>
#include <sstream>
#include <unordered_map>

#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"

#ifdef OMEGA_H_USE_GMSH
#include <esmfWrapper.h>
#endif  // OMEGA_H_USE_GMSH

namespace Omega_h {

namespace esmf {

namespace {

enum { //FIXME needed?
  GMSH_LINE = 1,
  GMSH_TRI = 2,
  GMSH_QUAD = 3,
  GMSH_TET = 4,
  GMSH_HEX = 5,
  GMSH_VERT = 15,
};

Int type_dim(Int type) { //FIXME needed?
  switch (type) {
    case GMSH_VERT:
      return 0;
    case GMSH_LINE:
      return 1;
    case GMSH_TRI:
    case GMSH_QUAD:
      return 2;
    case GMSH_TET:
    case GMSH_HEX:
      return 3;
  }
  Omega_h_fail(
      "omega_h can only accept linear simplices and hypercubes from Gmsh");
  OMEGA_H_NORETURN(-1);
}

Omega_h_Family type_family(Int type) {
  switch (type) {
    case GMSH_VERT:
    case GMSH_LINE:
    case GMSH_TRI:
    case GMSH_TET:
      return OMEGA_H_SIMPLEX;
    case GMSH_QUAD:
    case GMSH_HEX:
      return OMEGA_H_HYPERCUBE;
  }
  OMEGA_H_NORETURN(Omega_h_Family());
}

Int gmsh_type(Omega_h_Family family, Int dim) {
  switch (family) {
    case OMEGA_H_SIMPLEX:
      switch (dim) {
        case 0:
          return GMSH_VERT;
        case 1:
          return GMSH_LINE;
        case 2:
          return GMSH_TRI;
        case 3:
          return GMSH_TET;
      }
      return -1;
    case OMEGA_H_HYPERCUBE:
      switch (dim) {
        case 0:
          return GMSH_VERT;
        case 1:
          return GMSH_LINE;
        case 2:
          return GMSH_QUAD;
        case 3:
          return GMSH_HEX;
      }
      return -1;
    case OMEGA_H_MIXED:
      return -1;
  }
  return -1;
}

template <class T>
static void read(
    std::istream& stream, T& value, bool is_binary, bool needs_swapping) {
  if (is_binary) {
    binary::read_value(stream, value, needs_swapping);
  } else {
    stream >> value;
  }
}

}  // end anonymous namespace

Mesh read(std::istream& stream, CommPtr comm) {
  auto mesh = Mesh(comm->library());
  //determine how many parts the esmf mesh has and read on those ranks
  // - start with serial
  if (comm->rank() == 0) {
    esmfRead(stream, &mesh);
  }
  mesh.set_comm(comm);
  mesh.balance();
  return mesh;
}

Mesh read(filesystem::path const& filename, CommPtr comm) {
  std::ifstream file(filename.c_str());
  if (!file.is_open()) {
    Omega_h_fail("couldn't open \"%s\"\n", filename.c_str());
  }
  return esmf::read(file, comm);
}

}  // namespace esmf

}  // end namespace Omega_h
