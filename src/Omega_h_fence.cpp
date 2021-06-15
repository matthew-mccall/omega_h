#if defined(OMEGA_H_USE_SYCL)
#include <CL/sycl.hpp>
#include <dpct/dpct.hpp>
#endif
#include <Omega_h_fence.hpp>
#include <Omega_h_fail.hpp>

#ifdef OMEGA_H_USE_KOKKOS
#include <Omega_h_kokkos.hpp>
#endif

namespace Omega_h {

void fence() 
#if defined(OMEGA_H_USE_SYCL)
  try 
#endif
{
#if defined(OMEGA_H_USE_KOKKOS)
  Kokkos::fence();
#elif defined(OMEGA_H_USE_CUDA)
  auto const err = cudaDeviceSynchronize();
  OMEGA_H_CHECK(err == cudaSuccess);
#elif defined(OMEGA_H_USE_SYCL)
  dpct::get_current_device().queues_wait_and_throw();
#endif
}
#if defined(OMEGA_H_USE_SYCL)
catch (sycl::exception const &exc) {
  std::cerr << exc.what() << "Exception caught at file:" << __FILE__
            << ", line:" << __LINE__ << std::endl;
  std::exit(1);
}
#endif

}  // namespace Omega_h
