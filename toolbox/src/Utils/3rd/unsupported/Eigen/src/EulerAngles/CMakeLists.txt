# SPDX-FileCopyrightText: The Eigen Authors
# SPDX-License-Identifier: MPL-2.0

file(GLOB Eigen_EulerAngles_SRCS "*.h")

install(FILES
  ${Eigen_EulerAngles_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/unsupported/Eigen/src/EulerAngles COMPONENT Devel
  )
