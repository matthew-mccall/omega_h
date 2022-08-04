#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_vtk.hpp>

#ifdef OMEGA_H_USE_EGADS
#include <Omega_h_egads.hpp>
#endif


// Macro for checking errors in HIP API calls
#define hipErrorCheck(call)                                                                 \
do{                                                                                         \
    hipError_t hipErr = call;                                                               \
    if(hipSuccess != hipErr){                                                               \
        printf("HIP Error - %s:%d: '%s'\n", __FILE__, __LINE__, hipGetErrorString(hipErr)); \
        exit(0);                                                                            \
    }                                                                                       \
}while(0)


void printBinding(int rank) {
  // If ROCR_VISIBLE_DEVICES is set, capture visible GPUs
  const char* gpu_id_list; 
  const char* rocr_visible_devices = getenv("ROCR_VISIBLE_DEVICES");
  if(rocr_visible_devices == NULL){
    gpu_id_list = "N/A";
  }
  else{
    gpu_id_list = rocr_visible_devices;
  }

  std::string name = "crusher";

  // Find how many GPUs HIP runtime says are available
  int num_devices = 0;
  hipErrorCheck( hipGetDeviceCount(&num_devices) );

  int hwthread;
  int thread_id = 0;


  char busid[64];

  std::string busid_list = "";
  std::string rt_gpu_id_list = "";

  fprintf(stderr, "num_devices %d\n", num_devices);
  assert(num_devices == 1);
  // Loop over the GPUs available to each MPI rank
  for(int i=0; i<num_devices; i++){

    // Get the PCIBusId for each GPU and use it to query for UUID
    hipErrorCheck( hipDeviceGetPCIBusId(busid, 64, i) );

    // Concatenate per-MPIrank GPU info into strings for print
    if(i > 0) rt_gpu_id_list.append(",");
    rt_gpu_id_list.append(std::to_string(i));

    std::string temp_busid(busid);

    if(i > 0) busid_list.append(",");
    busid_list.append(temp_busid.substr(5,2));
  }

  thread_id = 0;
  hwthread = sched_getcpu();

  printf("MPI %03d - OMP %03d - HWT %03d - Node %s - RT_GPU_ID %s - GPU_ID %s - Bus_ID %s\n",
      rank, thread_id, hwthread, name.c_str(), rt_gpu_id_list.c_str(), gpu_id_list, busid_list.c_str());
}




int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  printBinding(lib.world()->rank());
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh_in.meshb");
  cmdline.add_arg<std::string>("mach_in.solb");
  cmdline.add_arg<std::string>("metric_in.solb");
#ifdef OMEGA_H_USE_EGADS
  cmdline.add_arg<std::string>("geom_in.egads");
#endif
  cmdline.add_arg<std::string>("out_prefix");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto mesh_path = cmdline.get<std::string>("mesh_in.meshb");
  auto mach_path = cmdline.get<std::string>("mach_in.solb");
  auto metric_path = cmdline.get<std::string>("metric_in.solb");
#ifdef OMEGA_H_USE_EGADS
  auto geom_path = cmdline.get<std::string>("geom_in.egads");
#endif
  auto out_prefix = cmdline.get<std::string>("out_prefix");
  Omega_h::Mesh mesh(&lib);
  Omega_h::meshb::read(&mesh, mesh_path.c_str());
  Omega_h::meshb::read_sol(&mesh, mach_path.c_str(), "mach");
  Omega_h::meshb::read_sol(&mesh, metric_path.c_str(), "original_metric");
#ifdef OMEGA_H_USE_EGADS
  auto geom = Omega_h::egads_load(geom_path);
#endif
  auto original_metrics =
      mesh.get_array<Omega_h::Real>(Omega_h::VERT, "original_metric");
  fprintf(stderr, "nverts %d nelms %d\n", mesh.nverts(), mesh.nelems());
  fprintf(stderr, "original_metrics size %d\n", original_metrics.size());
  auto graded_metrics =
      Omega_h::limit_metric_gradation(&mesh, original_metrics, 1.0);
  mesh.add_tag(Omega_h::VERT, "target_metric", Omega_h::symm_ncomps(mesh.dim()),
      graded_metrics);
  Omega_h::add_implied_metric_tag(&mesh);
  mesh.ask_qualities();
  auto opts = Omega_h::AdaptOpts(&mesh);
  opts.xfer_opts.type_map["mach"] = OMEGA_H_LINEAR_INTERP;
  opts.min_quality_allowed = 0.1;
#ifdef OMEGA_H_USE_EGADS
  opts.egads_model = geom;
#endif
  opts.verbosity = Omega_h::EXTRA_STATS;
  while (Omega_h::approach_metric(&mesh, opts)) {
    Omega_h::adapt(&mesh, opts);
  }
#ifdef OMEGA_H_USE_EGADS
  Omega_h::egads_free(geom);
#endif
  return 0;
}
