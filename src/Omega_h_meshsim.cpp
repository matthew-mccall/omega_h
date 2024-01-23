#include "Omega_h_file.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_mixedMesh.hpp"
#include "Omega_h_adj.hpp"

#include "SimModel.h"
#include "SimUtil.h"
#include "SimDiscrete.h"

namespace {
  int classId(pEntity e) {
    pGEntity g = EN_whatIn(e);
    assert(g);
    return GEN_tag(g);
  }

  int classType(pEntity e) {
    pGEntity g = EN_whatIn(e);
    assert(g);
    assert((0 <= GEN_type(g)) && (3 >= GEN_type(g)));
    return GEN_type(g);
  }

  int getNumber(pMeshNex nex, pEntity e, int order=0) {
    if(!MeshNex_hasNode(nex, e, order)) {
      fprintf(stderr, "ERROR: entity type %d id %d has no MeshNex node\n", EN_type(e), EN_id(e));
    }
    assert(MeshNex_hasNode(nex, e, order));
    return MeshNex_node(nex, e, order);
  }
}

namespace Omega_h {

namespace meshsim {

struct SimMeshInfo {
  int count_tet;
  int count_hex;
  int count_wedge;
  int count_pyramid;
  int count_tri;
  int count_quad;
  bool is_simplex;
  bool is_hypercube;
};

SimMeshInfo getSimMeshInfo(pMesh m) {
  SimMeshInfo info;

  RIter regions = M_regionIter(m);
  pRegion rgn;
  while ((rgn = (pRegion) RIter_next(regions))) {
    if (R_topoType(rgn) == Rtet) {
      info.count_tet += 1;
    }
    else if (R_topoType(rgn) == Rhex) {
      info.count_hex += 1;
    }
    else if (R_topoType(rgn) == Rwedge) {
      info.count_wedge += 1;
    }
    else if (R_topoType(rgn) == Rpyramid) {
      info.count_pyramid += 1;
    }
    else {
      Omega_h_fail("Region is not tet, hex, wedge, or pyramid \n");
    }
  }
  RIter_delete(regions);

  FIter faces = M_faceIter(m);
  pFace face;
  while ((face = (pFace) FIter_next(faces))) {
    if (F_numEdges(face) == 3) {
      info.count_tri += 1;
    }
    else if (F_numEdges(face) == 4) {
      info.count_quad += 1;
    }
    else {
      Omega_h_fail ("Face is neither tri nor quad \n");
    }
  }
  FIter_delete(faces);

  info.is_simplex = false;
  info.is_hypercube = false;
  if (info.count_hex == 0 && info.count_wedge == 0 && info.count_pyramid == 0) {
    if (info.count_tet == 0) {
      if (info.count_tri > 0) {
        info.is_simplex = true;
      }
    }
    else if (info.count_tet > 0) {
      info.is_simplex = true;
    }
    else {
      Omega_h_fail ("Invaild topology type\n");
    }
  }

  if (info.count_tet == 0 && info.count_wedge == 0 && info.count_pyramid == 0) {
    if (info.count_hex == 0) {
      if (info.count_quad > 0) {
        info.is_hypercube = true;
      }
    }
    else if (info.count_hex > 0) {
      info.is_hypercube = true;
    }
    else {
      Omega_h_fail ("Invaild topology type\n");
    }
  }

  fprintf(stderr, "tet=%d, hex=%d, wedge=%d, pyramid=%d\n",
         info.count_tet, info.count_hex, info.count_wedge, info.count_pyramid);
  fprintf(stderr, "tri=%d, quad=%d\n", info.count_tri, info.count_quad);
  return info;
}

struct SimMeshEntInfo {
  int maxDim;
  bool hasNumbering;
  std::vector<int> rgn_vertices[4];
  std::vector<int> ent_class_ids[4];
  std::vector<int> ent_class_dim[4];
  std::vector<int> ent_numbering;
  HostWrite<Real> host_coords;
  HostWrite<LO> host_class_ids_vtx;
  HostWrite<I8> host_class_dim_vtx;
  HostWrite<LO> host_class_ids_edge;
  HostWrite<I8> host_class_dim_edge;

  HostWrite<LO> host_e2v;

  SimMeshEntInfo(std::array<int,4> numEnts, bool hasNumbering_in) {
    hasNumbering = hasNumbering_in;
    maxDim = getMaxDim(numEnts);
    for(int dim=0; dim<4; dim++) {
      ent_class_ids[dim].reserve(numEnts[dim]);
      ent_class_dim[dim].reserve(numEnts[dim]);
    }
    if(hasNumbering) {
      ent_numbering.reserve(numEnts[0]);
    }
    host_coords = HostWrite<Real>(numEnts[0]*maxDim);

    host_class_ids_vtx = HostWrite<LO>(numEnts[0]);
    host_class_dim_vtx = HostWrite<I8>(numEnts[0]);

    host_class_ids_edge = HostWrite<LO>(numEnts[1]);
    host_class_dim_edge = HostWrite<I8>(numEnts[1]);

    host_e2v = HostWrite<LO>(numEnts[1]*2);
  }

  void readVerts(pMesh m,pMeshNex numbering) {
    VIter vertices = M_vertexIter(m);
    pVertex vtx;
    LO v = 0;
    while ((vtx = (pVertex) VIter_next(vertices))) {
      double xyz[3];
      V_coord(vtx,xyz);
      if( maxDim < 3 && xyz[2] != 0 )
        Omega_h_fail("The z coordinate must be zero for a 2d mesh!\n");
      for(int j=0; j<maxDim; j++) {
        host_coords[v * maxDim + j] = xyz[j];
      }
      ent_class_ids[0].push_back(classId(vtx));
      ent_class_dim[0].push_back(classType(vtx));
      if(hasNumbering) {
        ent_numbering.push_back(getNumber(numbering,vtx));
      }
      ++v;
    }
    VIter_delete(vertices);

    // WHY NOT STORE DIRECTORY INTO HostWrite???
    const int numVtx = M_numVertices(m);
    for (int i = 0; i < numVtx; ++i) {
      host_class_ids_vtx[i] = ent_class_ids[0][static_cast<std::size_t>(i)];
      host_class_dim_vtx[i] = ent_class_dim[0][static_cast<std::size_t>(i)];
    }
  }

  void readEdges(pMesh m) {
    std::vector<int> edge_vertices[1]; //TODO remove the '[1]'
    const int numEdges = M_numEdges(m);
    edge_vertices[0].reserve(numEdges*2);
    EIter edges = M_edgeIter(m);
    pEdge edge;
    int count_edge = 0;
    while ((edge = (pEdge) EIter_next(edges))) {
      double xyz[3];
      count_edge += 1;
      for(int j=0; j<2; ++j) {
        pVertex vtx = E_vertex(edge,j);
        edge_vertices[0].push_back(EN_id(vtx));
        V_coord(vtx,xyz);
      }
      ent_class_ids[1].push_back(classId(edge));
      ent_class_dim[1].push_back(classType(edge));
    }
    EIter_delete(edges);

    for (int i = 0; i < numEdges; ++i) {
      host_class_ids_edge[i] = ent_class_ids[1][static_cast<std::size_t>(i)];
      host_class_dim_edge[i] = ent_class_dim[1][static_cast<std::size_t>(i)];
    }

    for (Int i = 0; i < numEdges; ++i) {
      for (Int j = 0; j < 2; ++j) {
        host_e2v[i*2 + j] =
          edge_vertices[0][static_cast<std::size_t>(i*2 + j)];
      }
    }
  }

  struct MixedFaceClass {
    HostWrite<LO> triId;
    HostWrite<I8> triDim;
    std::vector<int> triVerts;
    HostWrite<LO> quadId;
    HostWrite<I8> quadDim;
    std::vector<int> quadVerts;
  };

  //mixed output:
  //tri2verts - face_vertices[0]
  //host_class_ids_tri
  //host_class_dim_tri
  //quads2verts - face_vertices[1]
  //host_class_ids_quad
  //host_class_dim_quad
  MixedFaceClass readMixedFaces(pMesh m, GO count_tri, GO count_quad) {
    std::vector<int> face_vertices[2];
    face_vertices[0].reserve(count_tri*3);
    face_vertices[1].reserve(count_quad*4);
    HostWrite<LO> host_class_ids_tri(count_tri);
    HostWrite<I8> host_class_dim_tri(count_tri);
    HostWrite<LO> host_class_ids_quad(count_quad);
    HostWrite<I8> host_class_dim_quad(count_quad);

    FIter faces = M_faceIter(m);
    pFace face;
    int triIdx = 0;
    int quadIdx = 0;
    while ((face = (pFace) FIter_next(faces))) {
      if (F_numEdges(face) == 3) {
        pVertex tri_vertex;
        pPList tri_vertices = F_vertices(face,1);
        assert (PList_size(tri_vertices) == 3);
        void *iter = 0;
        while ((tri_vertex = (pVertex) PList_next(tri_vertices, &iter))) {
          face_vertices[0].push_back(EN_id(tri_vertex));
        }
        PList_delete(tri_vertices);
        host_class_ids_tri[triIdx] = classId(face);
        host_class_dim_tri[triIdx] = classType(face);
        triIdx++;
      }
      else if (F_numEdges(face) == 4) {
        pVertex quad_vertex;
        pPList quad_vertices = F_vertices(face,1);
        assert (PList_size(quad_vertices) == 4);
        void *iter = 0;
        while ((quad_vertex = (pVertex) PList_next(quad_vertices, &iter))) {
          face_vertices[1].push_back(EN_id(quad_vertex));
        }
        PList_delete(quad_vertices);
        host_class_ids_quad[quadIdx] = classId(face);
        host_class_dim_quad[quadIdx] = classType(face);
        quadIdx++;
      }
      else {
        Omega_h_fail ("Face is neither tri nor quad \n");
      }
    }
    FIter_delete(faces);

    return MixedFaceClass(
             {host_class_ids_tri, host_class_dim_tri, face_vertices[0],
              host_class_ids_quad, host_class_dim_quad, face_vertices[1]}
           );
  }
 
  struct MonoFaceClass {
    HostWrite<LO> id;
    HostWrite<I8> dim;
    std::vector<int> verts;
  };

  //mono output:
  //==simplex==
  //tri2verts - face_vertices[0]
  //host_class_ids_face
  //host_class_dim_face
  //==hypecube==
  //quad2verts - face_vertices[1]
  //host_class_ids_face
  //host_class_dim_face
  MonoFaceClass readMonoTopoFaces(pMesh m, GO numFaces, LO vtxPerFace) {
    std::vector<int> face_vertices(numFaces*vtxPerFace);
    HostWrite<LO> host_face_class_ids(numFaces);
    HostWrite<I8> host_face_class_dim(numFaces);

    FIter faces = M_faceIter(m);
    pFace face;
    int faceIdx = 0;
    while ((face = (pFace) FIter_next(faces))) {
      assert(F_numEdges(face) == vtxPerFace);
      pVertex vertex;
      pPList vertices = F_vertices(face,1);
      assert (PList_size(vertices) == vtxPerFace);
      void *iter = 0;
      while ((vertex = (pVertex) PList_next(vertices, &iter))) {
        face_vertices.push_back(EN_id(vertex));
      }
      PList_delete(vertices);
      host_face_class_ids[faceIdx] = classId(face);
      host_face_class_dim[faceIdx] = classType(face);
      faceIdx++;
    }
    FIter_delete(faces);
    return MonoFaceClass({host_face_class_ids, host_face_class_dim, face_vertices});
  }

  private:
  SimMeshEntInfo();
  int getMaxDim(std::array<int,4> numEnts) {
    int max_dim;
    if (numEnts[3]) {
      max_dim = 3;
    } else if (numEnts[2]) {
      max_dim = 2;
    } else if (numEnts[1]) {
      max_dim = 1;
    } else {
      Omega_h_fail("There were no Elements of dimension higher than zero!\n");
    }
    return max_dim;
  }
};


void readMixed_internal(pMesh m, MixedMesh* mesh, SimMeshInfo info) {
  assert(!info.is_simplex && !info.is_hypercube);
  if(info.is_simplex || info.is_hypercube) {
      Omega_h_fail("Attempting to use the mixed mesh reader for a mono topology"
                   "mesh (family = simplex|hypercube)!\n");
  }

  const int numVtx = M_numVertices(m);
  const int numEdges = M_numEdges(m);
  const int numFaces = M_numFaces(m);
  const int numRegions = M_numRegions(m);

  const bool hasNumbering = false;

  SimMeshEntInfo simEnts({{numVtx,numEdges,numFaces,numRegions}}, hasNumbering);

  simEnts.readVerts(m,nullptr);

  mesh->set_dim(simEnts.maxDim);
  mesh->set_family(OMEGA_H_MIXED);

  mesh->set_verts_type(numVtx);
  mesh->add_coords_mix(simEnts.host_coords.write());
  mesh->add_tag<ClassId>(Topo_type::vertex, "class_id", 1,
                Read<ClassId>(simEnts.host_class_ids_vtx.write()));
  mesh->add_tag<I8>(Topo_type::vertex, "class_dim", 1,
                    Read<I8>(simEnts.host_class_dim_vtx.write()));

  simEnts.readEdges(m);

  auto ev2v = Read<LO>(simEnts.host_e2v.write());
  mesh->set_ents(Topo_type::edge, Topo_type::vertex, Adj(ev2v));
  mesh->template add_tag<ClassId>(Topo_type::edge, "class_id", 1,
                         Read<ClassId>(simEnts.host_class_ids_edge.write()));
  mesh->template add_tag<I8>(Topo_type::edge, "class_dim", 1,
                    Read<I8>(simEnts.host_class_dim_edge.write()));

  //process faces
  auto mixedFaceClass = simEnts.readMixedFaces(m, info.count_tri, info.count_quad);
  auto edge2vert = mesh->get_adj(Topo_type::edge, Topo_type::vertex);
  auto vert2edge = mesh->ask_up(Topo_type::vertex, Topo_type::edge);

  //// tris
  HostWrite<LO> host_tri2verts(info.count_tri*3);
  for (Int i = 0; i < info.count_tri; ++i) {
    for (Int j = 0; j < 3; ++j) {
      host_tri2verts[i*3 + j] =
          mixedFaceClass.triVerts[static_cast<std::size_t>(i*3 + j)];
    }
  }
  auto tri2verts = Read<LO>(host_tri2verts.write());
  auto down = reflect_down(tri2verts, edge2vert.ab2b, vert2edge,
      Topo_type::triangle, Topo_type::edge);
  mesh->set_ents(Topo_type::triangle, Topo_type::edge, down);
  mesh->template add_tag<ClassId>(Topo_type::triangle, "class_id", 1,
      Read<ClassId>(mixedFaceClass.triId.write()));
  mesh->template add_tag<I8>(Topo_type::triangle, "class_dim", 1,
      Read<I8>(mixedFaceClass.triDim.write()));

  //// quads
  HostWrite<LO> host_quad2verts(info.count_quad*4);
  for (Int i = 0; i < info.count_quad; ++i) {
    for (Int j = 0; j < 4; ++j) {
      host_quad2verts[i*4 + j] =
          mixedFaceClass.quadVerts[static_cast<std::size_t>(i*4 + j)];
    }
  }
  auto quad2verts = Read<LO>(host_quad2verts.write());
  down = reflect_down(quad2verts, edge2vert.ab2b, vert2edge,
      Topo_type::quadrilateral, Topo_type::edge);
  mesh->set_ents(Topo_type::quadrilateral, Topo_type::edge, down);
  mesh->template add_tag<ClassId>(Topo_type::quadrilateral, "class_id", 1,
      Read<ClassId>(mixedFaceClass.quadId.write()));
  mesh->template add_tag<I8>(Topo_type::quadrilateral, "class_dim", 1,
      Read<I8>(mixedFaceClass.quadDim.write()));

  //process regions - TODO
}


void read_internal(pMesh m, Mesh* mesh, pMeshNex numbering, SimMeshInfo info) {
  assert(info.is_simplex || info.is_hypercube);
  if(!info.is_simplex && !info.is_hypercube) {
      Omega_h_fail("Attempting to use the mono topology reader for a mixed"
                   "mesh (family = mixed)!\n");
  }

  const int numVtx = M_numVertices(m);
  const int numEdges = M_numEdges(m);
  const int numFaces = M_numFaces(m);
  const int numRegions = M_numRegions(m);

  const bool hasNumbering = false;

  SimMeshEntInfo simEnts({{numVtx,numEdges,numFaces,numRegions}}, hasNumbering);
  mesh->set_dim(simEnts.maxDim);

  //process verts
  simEnts.readVerts(m,nullptr);

  if (info.is_simplex) {
    mesh->set_family(OMEGA_H_SIMPLEX);
  } else if (info.is_hypercube){
    mesh->set_family(OMEGA_H_HYPERCUBE);
  }

  mesh->set_verts(numVtx);
  mesh->add_coords(simEnts.host_coords.write());
  mesh->add_tag<ClassId>(0, "class_id", 1,
      Read<ClassId>(simEnts.host_class_ids_vtx.write()));
  mesh->add_tag<I8>(0, "class_dim", 1,
      Read<I8>(simEnts.host_class_dim_vtx.write()));
  if(numbering) {
    HostWrite<LO> host_numbering_vtx(numVtx);
    for (int i = 0; i < numVtx; ++i)
      host_numbering_vtx[i] = simEnts.ent_numbering[static_cast<std::size_t>(i)];
    mesh->add_tag<LO>(0, "simNumbering", 1, Read<LO>(host_numbering_vtx.write()));
  }

  //process edges
  simEnts.readEdges(m);

  auto ev2v = Read<LO>(simEnts.host_e2v.write());
  mesh->set_ents(1, Adj(ev2v));
  mesh->add_tag<ClassId>(1, "class_id", 1,
                Read<ClassId>(simEnts.host_class_ids_edge.write()));
  mesh->add_tag<I8>(1, "class_dim", 1,
                    Read<I8>(simEnts.host_class_dim_edge.write()));

  //process faces
  if(info.is_simplex) {
    auto monoFaceClass = simEnts.readMonoTopoFaces(m, info.count_tri, 3);
    auto edge2vert = mesh->get_adj(1, 0);
    auto vert2edge = mesh->ask_up(0, 1);
    HostWrite<LO> host_tri2verts(info.count_tri*3);
    for (Int i = 0; i < info.count_tri; ++i) {
      for (Int j = 0; j < 3; ++j) {
        host_tri2verts[i*3 + j] =
          monoFaceClass.verts[static_cast<std::size_t>(i*3 + j)];
      }
    }
    auto tri2verts = Read<LO>(host_tri2verts.write());
    auto down = reflect_down(tri2verts, edge2vert.ab2b, vert2edge,
                             OMEGA_H_SIMPLEX, 2, 1);
    mesh->set_ents(2, down);
    mesh->add_tag<ClassId>(2, "class_id", 1,
                           Read<ClassId>(monoFaceClass.id.write()));
    mesh->add_tag<I8>(2, "class_dim", 1,
                      Read<I8>(monoFaceClass.dim.write()));
  }
  else {
    auto monoFaceClass = simEnts.readMonoTopoFaces(m, info.count_quad, 4);
    auto edge2vert = mesh->get_adj(1, 0);
    auto vert2edge = mesh->ask_up(0, 1);
    HostWrite<LO> host_quad2verts(info.count_quad*4);
    for (Int i = 0; i < info.count_quad; ++i) {
      for (Int j = 0; j < 4; ++j) {
        host_quad2verts[i*4 + j] =
          monoFaceClass.verts[static_cast<std::size_t>(i*4 + j)];
      }
    }
    auto quad2verts = Read<LO>(host_quad2verts.write());
    auto down = reflect_down(quad2verts, edge2vert.ab2b, vert2edge,
                        OMEGA_H_HYPERCUBE, 2, 1);
    mesh->set_ents(2, down);
    mesh->template add_tag<ClassId>(2, "class_id", 1,
                           Read<ClassId>(monoFaceClass.id.write()));
    mesh->template add_tag<I8>(2, "class_dim", 1,
                      Read<I8>(monoFaceClass.dim.write()));
  }

  //process regions - TODO
}
 
MixedMesh readMixedImpl(filesystem::path const& mesh_fname,
    filesystem::path const& mdl_fname,
    CommPtr comm) {
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimDiscrete_start(0);
  pNativeModel nm = NULL;
  pProgress p = NULL;
  pGModel g = GM_load(mdl_fname.c_str(), nm, p);
  pMesh m = M_load(mesh_fname.c_str(), g, p);
  auto simMeshInfo = getSimMeshInfo(m);
  auto mesh = MixedMesh(comm->library());
  mesh.set_comm(comm);
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  meshsim::readMixed_internal(m, &mesh, simMeshInfo);
  M_release(m);
  GM_release(g);
  SimDiscrete_stop(0);
  SimModel_stop();
  return mesh;
}

Mesh readImpl(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    filesystem::path const& numbering_fname, CommPtr comm) {
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimDiscrete_start(0);
  pNativeModel nm = NULL;
  pProgress p = NULL;
  pGModel g = GM_load(mdl_fname.c_str(), nm, p);
  pMesh m = M_load(mesh_fname.c_str(), g, p);
  auto simMeshInfo = getSimMeshInfo(m);
  const bool hasNumbering = (numbering_fname.native() != std::string(""));
  pMeshNex numbering = hasNumbering ? MeshNex_load(numbering_fname.c_str(), m) : nullptr;
  auto mesh = Mesh(comm->library());
  mesh.set_comm(comm);
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  meshsim::read_internal(m, &mesh, numbering, simMeshInfo);
  if(hasNumbering) MeshNex_delete(numbering);
  M_release(m);
  GM_release(g);
  SimDiscrete_stop(0);
  SimModel_stop();
  return mesh;
}

Mesh read(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    filesystem::path const& numbering_fname, CommPtr comm) {
  return readImpl(mesh_fname, mdl_fname, numbering_fname, comm);
}

Mesh read(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    CommPtr comm) {
  return readImpl(mesh_fname, mdl_fname, std::string(""), comm);
}

MixedMesh readMixed(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    CommPtr comm) {
  return readMixedImpl(mesh_fname, mdl_fname, comm);
}

}  // namespace meshsim

}  // end namespace Omega_h
