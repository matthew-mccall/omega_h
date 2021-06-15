#if defined(OMEGA_H_USE_SYCL)
#include <CL/sycl.hpp>
#include <dpct/dpct.hpp>
#endif
#include <Omega_h_fail.hpp>
#include <Omega_h_malloc.hpp>
#include <Omega_h_pool.hpp>
#include <Omega_h_profile.hpp>
#include <cstdlib>

namespace Omega_h {

void *device_malloc(std::size_t size)
#if defined(OMEGA_H_USE_SYCL)
  try 
#endif
{
  OMEGA_H_TIME_FUNCTION;
#ifdef OMEGA_H_USE_CUDA
  void* tmp_ptr;
  auto cuda_malloc_size = size;
  if (cuda_malloc_size < 1) cuda_malloc_size = 1;
  auto const err = cudaMalloc(&tmp_ptr, cuda_malloc_size);
  if (err == cudaErrorMemoryAllocation) return nullptr;
  OMEGA_H_CHECK(err == cudaSuccess);
  return tmp_ptr;
#elif defined(OMEGA_H_USE_SYCL)
  void* tmp_ptr;
  auto sycl_malloc_size = size;
  if (sycl_malloc_size < 1) sycl_malloc_size = 1;
  tmp_ptr = (void *)sycl::malloc_device(
                        sycl_malloc_size, dpct::get_default_queue());
  return tmp_ptr;
#else
  return ::std::malloc(size);
#endif
}
#if defined(OMEGA_H_USE_SYCL)
catch (sycl::exception const &exc) {
  std::cerr << exc.what() << "Exception caught at file:" << __FILE__
            << ", line:" << __LINE__ << std::endl;
  std::exit(1);
}
#endif

void device_free(void *ptr, std::size_t)
#if defined(OMEGA_H_USE_SYCL)
  try 
#endif
{
  OMEGA_H_TIME_FUNCTION;
#ifdef OMEGA_H_USE_CUDA
  auto const err = cudaFree(ptr);
  OMEGA_H_CHECK(err == cudaSuccess);
#elif defined(OMEGA_H_USE_SYCL)
  sycl::free(ptr, dpct::get_default_queue());
#else
  ::std::free(ptr);
#endif
}
#if defined(OMEGA_H_USE_SYCL)
catch (sycl::exception const &exc) {
  std::cerr << exc.what() << "Exception caught at file:" << __FILE__
            << ", line:" << __LINE__ << std::endl;
  std::exit(1);
}
#endif

void* host_malloc(std::size_t size)
#if defined(OMEGA_H_USE_SYCL)
  try 
#endif
{
  OMEGA_H_TIME_FUNCTION;
#ifdef OMEGA_H_USE_CUDA
  void* tmp_ptr;
  auto cuda_malloc_size = size;
  if (cuda_malloc_size < 1) cuda_malloc_size = 1;
  auto const err = cudaMallocHost(&tmp_ptr, cuda_malloc_size);
  if (err == cudaErrorMemoryAllocation) return nullptr;
  OMEGA_H_CHECK(err == cudaSuccess);
  return tmp_ptr;
#elif defined(OMEGA_H_USE_SYCL)
  void* tmp_ptr;
  auto sycl_malloc_size = size;
  if (sycl_malloc_size < 1) sycl_malloc_size = 1;
  tmp_ptr = (void *)sycl::malloc_host(
                        sycl_malloc_size, dpct::get_default_queue());
  return tmp_ptr;
#else
  return ::std::malloc(size);
#endif
}
#if defined(OMEGA_H_USE_SYCL)
catch (sycl::exception const &exc) {
  std::cerr << exc.what() << "Exception caught at file:" << __FILE__
            << ", line:" << __LINE__ << std::endl;
  std::exit(1);
}
#endif

void host_free(void *ptr, std::size_t)
#if defined(OMEGA_H_USE_SYCL)
  try
#endif
{
  OMEGA_H_TIME_FUNCTION;
#ifdef OMEGA_H_USE_CUDA
  auto const err = cudaFreeHost(ptr);
  OMEGA_H_CHECK(err == cudaSuccess);
#elif defined(OMEGA_H_USE_SYCL)
  sycl::free(ptr, dpct::get_default_queue());
#else
  ::std::free(ptr);
#endif
}
#if defined(OMEGA_H_USE_SYCL)
catch (sycl::exception const &exc) {
  std::cerr << exc.what() << "Exception caught at file:" << __FILE__
            << ", line:" << __LINE__ << std::endl;
  std::exit(1);
}
#endif

static Pool* device_pool = nullptr;
static Pool* host_pool = nullptr;

void enable_pooling() {
  device_pool = new Pool(device_malloc, device_free);
  host_pool = new Pool(host_malloc, host_free);
}

void disable_pooling() {
  delete device_pool;
  delete host_pool;
  device_pool = nullptr;
  host_pool = nullptr;
}

void* maybe_pooled_device_malloc(std::size_t size) {
  if (device_pool) return allocate(*device_pool, size);
  return device_malloc(size);
}

void maybe_pooled_device_free(void* ptr, std::size_t size) {
  if (device_pool)
    deallocate(*device_pool, ptr, size);
  else
    device_free(ptr, size);
}

void* maybe_pooled_host_malloc(std::size_t size) {
  if (host_pool) return allocate(*host_pool, size);
  return host_malloc(size);
}

void maybe_pooled_host_free(void* ptr, std::size_t size) {
  if (host_pool)
    deallocate(*host_pool, ptr, size);
  else
    host_free(ptr, size);
}
}  // namespace Omega_h
