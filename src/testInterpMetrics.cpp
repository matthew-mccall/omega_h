#include <Omega_h_array_ops.hpp>
#include <Omega_h_eigen.hpp>
#include <Omega_h_lie.hpp>
#include <Omega_h_metric_intersect.hpp>
#include <Omega_h_most_normal.hpp>
#include <Omega_h_shape.hpp>
#include <Omega_h_svd.hpp>

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(std::string(lib.version()) == OMEGA_H_SEMVER);
  {
    auto a = repeat_symm(
        4, compose_metric(identity_matrix<2, 2>(), vector_2(1.0 / 100.0, 1.0)));
    auto b = repeat_symm(
        4, compose_metric(identity_matrix<2, 2>(), vector_2(1.0, 1.0)));
    auto c = interpolate_between_metrics(4, a, b, 0.0);
    OMEGA_H_CHECK(are_close(a, c));
    c = interpolate_between_metrics(4, a, b, 1.0);
    OMEGA_H_CHECK(are_close(b, c));
  }
  return 0;
}
