#include <Omega_h_adapt.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_profile.hpp>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("input.osh");
  cmdline.add_arg<std::string>("output.osh");
  cmdline.add_arg<std::string>("output.vtk");
  auto const world = lib.world();
  if (!cmdline.parse_final(world, &argc, argv)) {
    return -1;
  }
  Omega_h::ScopedTimer scoped_timer("main");
  auto const world_size = world->size();
  auto const world_rank = world->rank();
  auto const inpath = cmdline.get<std::string>("input.osh");
  auto const outpath = cmdline.get<std::string>("output.osh");
  auto const outvtkpath = cmdline.get<std::string>("output.vtk");
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[1], lib.world(), &mesh);
  Omega_h::AdaptOpts opts(&mesh);
  auto metrics = mesh.get_array<double>(0, "size");
  auto const metric_ncomps =
    Omega_h::divide_no_remainder(metrics.size(), mesh.nverts());
  mesh.add_tag(0, "metric", metric_ncomps, metrics);
  if (world_rank == 0) std::cout << "adapting to metric\n";
  Omega_h::adapt(&mesh, opts);
  auto nelems = mesh.nglobal_ents(mesh.dim());
  if (world_rank == 0)
    std::cout << "mesh now has " << nelems << " total elements\n";
  Omega_h::vtk::write_parallel(outvtkpath, &mesh, mesh.dim());
  return 0;
}
