#include "Omega_h_swap3d.hpp"

#include "Omega_h_for.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_swap3d_loop.hpp"
#include "Omega_h_swap3d_tables.hpp"

namespace Omega_h {


HostFew<LOs, 4> swap3d_topology(Mesh* mesh, LOs keys2edges,
    Read<I8> edge_configs, HostFew<LOs, 4> keys2prods) {
  auto edges2tets = mesh->ask_up(EDGE, REGION);
  auto edges2edge_tets = edges2tets.a2ab;
  auto edge_tets2tets = edges2tets.ab2b;
  auto edge_tet_codes = edges2tets.codes;
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto tet_verts2verts = mesh->ask_verts_of(REGION);
  HostFew<Write<LO>, 4> prod_verts2verts_w;
  for (Int prod_dim = EDGE; prod_dim <= REGION; ++prod_dim) {
    prod_verts2verts_w[prod_dim] =
        Write<LO>(keys2prods[prod_dim].last() * Int(prod_dim + 1));
  }
  auto nkeys = keys2edges.size();
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto edge = keys2edges[key];
    auto config = edge_configs[edge];
    auto loop = swap3d::find_loop(edges2edge_tets, edge_tets2tets,
        edge_tet_codes, edge_verts2verts, tet_verts2verts, edge);
    auto nplane_tris = swap3d::swap_mesh_sizes[loop.size];
    auto nplane_edges = swap3d::nedges[loop.size];
    for (Int plane_edge = 0; plane_edge < nplane_edges; ++plane_edge) {
      auto unique_edge = swap3d::edges2unique[loop.size][config][plane_edge];
      Few<LO, 2> plane_edge_verts;
      for (Int pev = 0; pev < 2; ++pev) {
        auto loop_vert = swap3d::unique_edges[loop.size][unique_edge][pev];
        auto vert = loop.loop_verts2verts[loop_vert];
        plane_edge_verts[pev] = vert;
      }
      auto prod_edge = keys2prods[EDGE][key] + plane_edge;
      for (Int pev = 0; pev < 2; ++pev) {
        prod_verts2verts_w[EDGE][prod_edge * 2 + pev] = plane_edge_verts[pev];
      }
    }
  };
  parallel_for(1, f, "swap3d_topology");
  HostFew<LOs, 4> prod_verts2verts;
  return prod_verts2verts;
}

}  // end namespace Omega_h
