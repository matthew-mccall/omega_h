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
  Omega_h::binary::read(inpath, lib.world(), &mesh);

  mesh.ask_sizes();
  Omega_h::vtk::write_parallel("size.vtk", &mesh, mesh.dim());

  //enable ghosting for adaptation
  mesh.set_parting(OMEGA_H_GHOSTED);

  //apapt options
  Omega_h::AdaptOpts opts(&mesh);
  opts.verbosity = Omega_h::EXTRA_STATS;

  //used to define 'target_metric' needed by approach_metric
  auto metrics = mesh.get_array<double>(0, "size");
  mesh.remove_tag(0,"size"); //"size" appears to be a reserved keyword  ... at least for elements

  //adapt - this creates the needed 'metric' and 'target_metric' tags
  Omega_h::grade_fix_adapt(&mesh,opts,metrics,/*verbose*/true);
  
  Omega_h::vtk::write_parallel(outvtkpath, &mesh, mesh.dim());
  auto nelems = mesh.nglobal_ents(mesh.dim());

  if (world_rank == 0)
    std::cout << "mesh now has " << nelems << " total elements\n";
  return 0;
}
