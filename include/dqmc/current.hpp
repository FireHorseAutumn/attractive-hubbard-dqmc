#pragma once

#include "dqmc/dqmc.hpp"

namespace dqmc {

CurrentResult exact_tiny_current_response(const Matrix& kinetic, double abs_u,
                                          double beta, double t, int lx, int ly,
                                          int qx, int qy);

namespace test_support {

Matrix inverse_identity_plus(const Matrix& matrix);
Matrix inverse_inverse_plus(const Matrix& u, const Matrix& v);
Matrix product_ldr(const std::vector<Matrix>& factors);

}  // namespace test_support

}  // namespace dqmc
