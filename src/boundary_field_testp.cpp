#include <iostream>
#include <string>

#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_array_ops.hpp"
#include <Omega_h_adapt.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_timer.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_build.hpp>

using namespace Omega_h;

template <Int dim>
static void set_target_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto target_metrics_w = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto z = coords[v * dim + (dim - 1)];
    auto h = Vector<dim>();
    for (Int i = 0; i < dim - 1; ++i) h[i] = 0.1;
    h[dim - 1] = 0.001 + 0.198 * std::abs(z - 0.5);
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  parallel_for(mesh->nverts(), f);
  mesh->set_tag(VERT, "target_metric", Reals(target_metrics_w));
}

template <Int dim>
void run_case(Mesh* mesh, char const* vtk_path) {
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = get_implied_metrics(mesh);
  mesh->add_tag(0, "metric", symm_ncomps(dim), implied_metrics);
  mesh->add_tag<Real>(VERT, "target_metric", symm_ncomps(dim));
  set_target_metric<dim>(mesh);

  OMEGA_H_CHECK(mesh->has_tag(0, "metric"));

  mesh->set_parting(OMEGA_H_ELEM_BASED);

  OMEGA_H_CHECK(mesh->has_tag(0, "metric"));

  mesh->ask_lengths();
  mesh->ask_qualities();
  vtk::FullWriter writer;
  if (vtk_path) {
    writer = vtk::FullWriter(vtk_path, mesh);
    writer.write();
  }

  auto opts = AdaptOpts(mesh);
  opts.verbosity = EXTRA_STATS;
  opts.length_histogram_max = 2.0;
  opts.max_length_allowed = opts.max_length_desired * 2.0;
  add_boundaryField_transferMap(&opts, "field1", OMEGA_H_LINEAR_INTERP);
  //add_boundaryField_transferMap(&opts, "field2", OMEGA_H_LINEAR_INTERP);
  add_boundaryField_transferMap(&opts, "field2", OMEGA_H_INHERIT);
  add_boundaryField_transferMap(&opts, "field3", OMEGA_H_INHERIT);
  //opts.xfer_opts.type_map["magnetic_face_flux"] = OMEGA_H_POINTWISE;
  //opts.xfer_opts.integral_map["magnetic_face_flux"] = "mass";
  
  //add_boundaryField_transferMap(&opts, "field3", OMEGA_H_INHERIT);
  //add_boundaryField_transferMap(&opts, "field3", OMEGA_H_CONSERVE);
  //add_boundaryField_integralMap(&opts, "field3", "mass");
  add_boundaryField_transferMap(&opts, "field4", OMEGA_H_POINTWISE);
  Now t0 = now();
  while (approach_metric(mesh, opts)) {
    adapt(mesh, opts);
    if (mesh->has_tag(VERT, "target_metric")) set_target_metric<dim>(mesh);

    auto face_rc = mesh->ask_revClass(2);
    auto face_a2abSize = face_rc.a2ab.size();
    OMEGA_H_CHECK(face_a2abSize);
    auto ab2b = face_rc.ab2b;
    auto a2ab = face_rc.a2ab;
    auto nbface = ab2b.size();
    Write<Real> valsf_new(nbface, 12);

    auto f_new = OMEGA_H_LAMBDA(LO gf) {
      if ((gf == 16) || (gf == 22)) {
        auto start = a2ab[gf];
        auto end = a2ab[gf+1];
        for (int index = start; index < end; ++index) {
          valsf_new[index] = gf;
        }
      }
    };
    parallel_for(face_a2abSize-1, f_new);

    mesh->set_boundaryField_array<Real>(2, "field3", Read<Real>(valsf_new));
    if (vtk_path) writer.write();
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";

}

void test_3d(Library *lib) {

  auto mesh = Mesh(lib);
  binary::read ("./../../omega_h/meshes/box_3d_2p.osh",
                lib->world(), &mesh);

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

  mesh.add_boundaryField<Real>(1, "field2", 1);
  const auto rank = lib->world()->rank();
  if ((!rank) || (rank == 1)) {
    Write<Real> vals(nbedge, 100);
    Read<Real> vals_r(vals);
    mesh.set_boundaryField_array<Real>(1, "field2", vals_r);
  }
  else {
    Write<Real> vals(nbedge, 50.45632);
    Read<Real> vals_r(vals);
    mesh.set_boundaryField_array<Real>(1, "field2", vals_r);
  }

  mesh.add_boundaryField<Real>(2, "field3", 1);
  Write<Real> valsf(nbface, 12);

  auto ab2b = face_rc.ab2b;
  auto a2ab = face_rc.a2ab;
  auto f = OMEGA_H_LAMBDA(LO gf) {
    if ((gf == 16) || (gf == 22)) {
      auto start = a2ab[gf];
      auto end = a2ab[gf+1];
      for (int index = start; index < end; ++index) {
        valsf[index] = gf;
      }
    }
  };
  parallel_for(face_a2abSize-1, f);

  Read<Real> valsf_r(valsf);
  mesh.set_boundaryField_array<Real>(2, "field3", valsf_r);
  Write<Real> vals_allFace(mesh.nfaces(), 12.5);
  mesh.add_tag<Real>(2, "magnetic_face_flux", 1, Read<Real>(vals_allFace));
 
  mesh.add_boundaryField<Real>(3, "field4", 1);
  Write<Real> valsr(nbreg, 100);
  Read<Real> valsr_r(valsr);
  mesh.set_boundaryField_array<Real>(3, "field4", valsr_r);

  vtk::write_parallel ("./../../omega_h/meshes/box_3d_2p.vtk",
                   &mesh);

  auto new_mesh = Mesh(lib);
  binary::read ("./../../omega_h/meshes/box_3d_2p.osh",
                lib->world(), &new_mesh);
  auto new_bField = new_mesh.get_boundaryField_array<Real>(0, "field1");
  auto vals_r = mesh.get_boundaryField_array<Real>(0, "field1");
  OMEGA_H_CHECK(new_bField == vals_r);
  auto nverts = mesh.nverts();
  OMEGA_H_CHECK(new_bField.size() <= nverts);

  mesh.sync_boundaryField(1, "field2");
  vtk::write_parallel
    ("./../../omega_h/meshes/box_3d_2p_sync.vtk", &mesh);
  mesh.reduce_boundaryField(1, "field2", OMEGA_H_SUM);
  vtk::write_parallel
    ("./../../omega_h/meshes/box_3d_2p_reduce.vtk", &mesh);

  vtk::write_parallel ("./../../omega_h/meshes/box_3d_foo_2p.vtk",
                   &mesh);
  run_case<3>(&mesh, "./../../omega_h/meshes/adapt/box3d_2p.vtk");

  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

  test_3d(&lib);

  return 0;
}
