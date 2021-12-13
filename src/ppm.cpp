#include "Omega_h_array_ops.hpp"
#include "Omega_h_compare.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace Omega_h;

#ifdef OMEGA_H_USE_ZLIB
#endif

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(std::string(lib.version()) == OMEGA_H_SEMVER);
  if (lib.world()->size() == 1) {
    /*
    std::ifstream in_str(ppm_file);
    if (!in_str) {
      fprintf(stderr, "[ERROR] Cannot open file %s\n", ppm_file);
      return;
    }

#ifdef OMEGA_H_USE_ZLIB
    bool compress = true;
#else
    bool compress = false;
#endif
    bool swap = !is_little_endian_cpu();
    Omega_h::I8 version;
    Omega_h::binary::read_value(in_str, version, swap);
    Omega_h::I8 is_full_mesh_int;
    Omega_h::binary::read_value(in_str, is_full_mesh_int, swap);
    int num_entities[4];
    int num_cores[4];
    LOs buffered_parts[4];
    LOs offset_ents_per_rank_per_dim[4];
    LOs ent_to_comm_arr_index_per_dim[4];
    LOs is_complete_part[4];

    for (int i = 0; i < 4; ++i) {
      if (version >= 2) {
        //Read num_entites
        Omega_h::binary::read_value(in_str, num_entities[i], swap);
      }
      //Read num_cores
      Omega_h::binary::read_value(in_str, num_cores[i], swap);
      //Read buffered_parts
      Omega_h::binary::read_array(in_str, buffered_parts[i], compress, swap);
      //Read offset_ents_per_rank_per_dim
      Omega_h::binary::read_array(in_str, offset_ents_per_rank_per_dim[i],
                                   compress, swap);
      //Read ent_to_comm_arr_index_per_dim
      Omega_h::binary::read_array(in_str, ent_to_comm_arr_index_per_dim[i],
                                   compress, swap);
      //Read is complete part
      assert(cudaSuccess == cudaDeviceSynchronize());
      fprintf(stderr, "0.1\n");
      Omega_h::binary::read_array(in_str, is_complete_part[i], compress, swap);
      fprintf(stderr, "0.2\n");
      ////read num_bounds
      //Omega_h::binary::read_value(in_str, num_bounds[i], swap);
      ////read num_boundaries
      //Omega_h::binary::read_value(in_str, num_boundaries[i], swap);
      ////read boundary_parts
      //Omega_h::binary::read_array(in_str, boundary_parts[i], compress, swap);
      ////read offset_bounded_per_dim
      //Omega_h::binary::read_array(in_str, offset_bounded_per_dim[i],
      //                            compress, swap);
      ////read bounded_ent_ids
      //Omega_h::binary::read_array(in_str, bounded_ent_ids[i], compress, swap);
    }
*/

  }
  return 0;
}
