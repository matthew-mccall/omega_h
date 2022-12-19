#include <Omega_h_library.hpp>

extern "C" void esmfInit();
extern "C" void esmfFinalize();
extern "C" void esmfTestMeshGet();

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  esmfInit();
  esmfTestMeshGet();
  esmfFinalize();
  return 0;
}
