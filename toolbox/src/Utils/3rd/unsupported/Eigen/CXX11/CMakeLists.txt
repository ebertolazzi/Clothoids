# Forwarding headers for backward compatibility.
# SPDX-FileCopyrightText: The Eigen Authors
# SPDX-License-Identifier: MPL-2.0
set(Eigen_CXX11_HEADERS Tensor TensorSymmetry ThreadPool)

install(FILES
  ${Eigen_CXX11_HEADERS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/unsupported/Eigen/CXX11 COMPONENT Devel
  )
