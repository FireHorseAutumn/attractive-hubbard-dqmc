#pragma once

#include "dqmc/dqmc.hpp"

#include <Eigen/LU>
#include <Eigen/QR>

namespace dqmc::detail {

struct Ldr {
  Matrix left;
  Vector diagonal;
  Matrix right;
};

inline Ldr identity_ldr(int n) {
  return {Matrix::Identity(n, n), Vector::Ones(n), Matrix::Identity(n, n)};
}

inline Ldr make_ldr(const Matrix& matrix) {
  const Eigen::FullPivHouseholderQR<Matrix> qr(matrix);
  Matrix r = qr.matrixQR().triangularView<Eigen::Upper>();
  const Vector diagonal = r.diagonal().cwiseAbs();
  r = diagonal.cwiseInverse().asDiagonal() * r;
  r = r * qr.colsPermutation().inverse();
  return {qr.matrixQ(), diagonal, r};
}

inline Ldr prepend_ldr(const Ldr& previous, const Matrix& b) {
  const Matrix w = b * previous.left * previous.diagonal.asDiagonal();
  const Eigen::HouseholderQR<Matrix> qr(w);
  Matrix r = qr.matrixQR().triangularView<Eigen::Upper>();
  const Vector diagonal = r.diagonal().cwiseAbs();
  r = diagonal.cwiseInverse().asDiagonal() * r * previous.right;
  return {qr.householderQ(), diagonal, r};
}

inline Matrix pivoted_inverse(const Matrix& matrix) {
  return Eigen::FullPivLU<Matrix>(matrix).inverse();
}

inline Matrix inverse_inverse_plus(const Ldr& u, const Ldr& v) {
  const Vector du_min = u.diagonal.cwiseMin(1.0);
  const Vector du_max = u.diagonal.cwiseMax(1.0);
  const Vector dv_min = v.diagonal.cwiseMin(1.0);
  const Vector dv_max = v.diagonal.cwiseMax(1.0);
  const Matrix rv_dv_max =
      pivoted_inverse(v.right) * dv_max.cwiseInverse().asDiagonal();
  const Matrix middle = du_max.cwiseInverse().asDiagonal() * u.left.transpose() *
                          rv_dv_max +
                        du_min.asDiagonal() * u.right * v.left * dv_min.asDiagonal();
  return rv_dv_max * pivoted_inverse(middle) * du_min.asDiagonal() * u.right;
}

inline Matrix inverse_identity_plus(const Ldr& a) {
  const Vector d_min = a.diagonal.cwiseMin(1.0);
  const Vector d_max = a.diagonal.cwiseMax(1.0);
  const Matrix r_d_max =
      pivoted_inverse(a.right) * d_max.cwiseInverse().asDiagonal();
  const Matrix middle = r_d_max + a.left * d_min.asDiagonal();
  return r_d_max * pivoted_inverse(middle);
}

}  // namespace dqmc::detail
