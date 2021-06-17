#ifndef OMEGA_H_ATOMICS_HPP
#define OMEGA_H_ATOMICS_HPP

#if defined(OMEGA_H_USE_SYCL)
#endif
#include <Omega_h_config.h>

#if defined(OMEGA_H_USE_KOKKOS)
#include <Kokkos_Core.hpp>
#endif

namespace Omega_h {

OMEGA_H_DEVICE int atomic_fetch_add(int* const dest, const int val) {
#if defined(OMEGA_H_USE_KOKKOS)
  return Kokkos::atomic_fetch_add(dest, val);
#elif defined(OMEGA_H_USE_OPENMP)
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"
#endif
  int oldval;
#pragma omp atomic capture
  {
    oldval = *dest;
    *dest += val;
  }
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
  return oldval;
#elif defined(OMEGA_H_USE_CUDA)
  return atomicAdd(dest, val);
#elif defined(OMEGA_H_USE_SYCL)
  /*
  DPCT1039:0: The generated code assumes that "dest" points to the global memory
  address space. If it points to a local memory address space, replace
  "dpct::atomic_fetch_add" with "dpct::atomic_fetch_add<int,
  sycl::access::address_space::local_space>".
  */
  return sycl::atomic<int>(sycl::global_ptr<int>(dest)).fetch_add(val);
#else
  int oldval = *dest;
  *dest += val;
  return oldval;
#endif
}

OMEGA_H_DEVICE void atomic_increment(int* const dest) {
#if defined(OMEGA_H_USE_OPENMP) || defined(OMEGA_H_USE_CUDA)
  atomic_fetch_add(dest, 1);
#else
  ++(*dest);
#endif
}

OMEGA_H_DEVICE void atomic_add(int* const dest, const int val) {
#if defined(OMEGA_H_USE_OPENMP) || defined(OMEGA_H_USE_CUDA)
  atomic_fetch_add(dest, val);
#else
  *dest += val;
#endif
}

}  // end namespace Omega_h

#endif
