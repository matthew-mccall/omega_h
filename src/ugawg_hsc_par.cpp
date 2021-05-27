#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_vtk.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh_in.osh");
  cmdline.add_arg<std::string>("out_prefix");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto mesh_path = cmdline.get<std::string>("mesh_in.osh");
  auto out_prefix = cmdline.get<std::string>("out_prefix");

  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(mesh_path, world, &mesh);
  world->barrier();

  auto imb = mesh.imbalance();
  if (!mesh.comm()->rank())
    fprintf(stderr, "mesh element imbalance %.2f\n", imb);

  auto nElms = mesh.nelems();
  auto rank = world->rank();
  fprintf(stderr, "%d elements %d\n", rank, nElms);

  auto target_metric =
      mesh.get_array<Omega_h::Real>(Omega_h::VERT, "target_metric");
  Omega_h::add_implied_metric_tag(&mesh);
  mesh.ask_qualities();
  auto opts = Omega_h::AdaptOpts(&mesh);
  opts.xfer_opts.type_map["mach"] = OMEGA_H_LINEAR_INTERP;
  opts.min_quality_allowed = 0.1;

  auto vtx_rc = mesh.ask_revClass(0);
  auto vert_boundary_ids = (mesh.ask_revClass(0)).ab2b;
  auto nbvert = vert_boundary_ids.size();
  OMEGA_H_CHECK (nbvert < mesh.nverts());
  auto edge_boundary_ids = (mesh.ask_revClass(1)).ab2b;
  auto nbedge = edge_boundary_ids.size();
  auto face_rc = mesh.ask_revClass(2);
  auto face_a2abSize = face_rc.a2ab.size();
  OMEGA_H_CHECK(face_a2abSize);
  auto face_boundary_ids = (mesh.ask_revClass(2)).ab2b;
  auto nbface = face_boundary_ids.size();
  auto reg_boundary_ids = (mesh.ask_revClass(3)).ab2b;
  auto nbreg = reg_boundary_ids.size();

  mesh.add_boundaryField<Omega_h::LO>(0, "field", 1);
  Omega_h::Write<Omega_h::LO> valsr_v(nbvert, 100);
  mesh.set_boundaryField_array(0, "field", Omega_h::Read<Omega_h::LO>(valsr_v));
  Omega_h::Write<Omega_h::Real> vals_real(nbvert, 50.2523232);
  mesh.add_boundaryField<Omega_h::Real>(0, "field1", 1, Omega_h::Read<Omega_h::Real>(vals_real));
  mesh.add_boundaryField<Omega_h::LO>(1, "field", 1);
  Omega_h::Write<Omega_h::LO> vals(nbedge, 100);
  Omega_h::Read<Omega_h::LO> vals_r(vals);
  mesh.set_boundaryField_array(1, "field", vals_r);
  mesh.add_boundaryField<Omega_h::LO>(2, "field", 1);
  Omega_h::Write<Omega_h::LO> valsf(nbface, 12);
  mesh.set_boundaryField_array(2, "field", Omega_h::Read<Omega_h::LO>(valsf));
  mesh.add_boundaryField<Omega_h::LO>(3, "field", 1);
  Omega_h::Write<Omega_h::LO> valsr(nbreg, 100);
  Omega_h::Read<Omega_h::LO> valsr_r(valsr);
  mesh.set_boundaryField_array(3, "field", valsr_r);

  add_boundaryField_transferMap(&opts, "field", OMEGA_H_INHERIT);
  add_boundaryField_transferMap(&opts, "field1", OMEGA_H_LINEAR_INTERP);

  opts.verbosity = Omega_h::EXTRA_STATS;
  while (Omega_h::approach_metric(&mesh, opts)) {
    Omega_h::adapt(&mesh, opts);
  }
  return 0;
}
