#ifndef OMEGA_H_KOKKOS_HPP
#define OMEGA_H_KOKKOS_HPP

#include <Omega_h_config.h>

#include <Omega_h_macros.h>

OMEGA_H_SYSTEM_HEADER

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include <Kokkos_Core.hpp>

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#if defined(KOKKOS_HAVE_CUDA) && !defined(OMEGA_H_USE_CUDA)
#error "Kokkos has CUDA, please reconfigure with Omega_h_USE_CUDA=ON"
#endif

namespace Omega_h {
using ExecSpace = Kokkos::DefaultExecutionSpace;
using StaticSched = Kokkos::Schedule<Kokkos::Static>;
using Policy = Kokkos::RangePolicy<ExecSpace, StaticSched, Omega_h::LO>;

inline Policy policy(LO n) { return Policy(0, n); }
}

#endif
