// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2014 Benoit Steiner <benoit.steiner.goog@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#include "main.h"

#include <Eigen/Tensor>

using Eigen::Tensor;

void test_simple_patch() {
  Tensor<float, 4> tensor(2, 3, 5, 7);
  tensor.setRandom();
  Tensor<float, 4, RowMajor> tensor_row_major = tensor.swap_layout();
  VERIFY_IS_EQUAL(tensor.dimension(0), tensor_row_major.dimension(3));
  VERIFY_IS_EQUAL(tensor.dimension(1), tensor_row_major.dimension(2));
  VERIFY_IS_EQUAL(tensor.dimension(2), tensor_row_major.dimension(1));
  VERIFY_IS_EQUAL(tensor.dimension(3), tensor_row_major.dimension(0));

  // Single pixel patch: ColMajor
  Tensor<float, 5> single_pixel_patch;
  single_pixel_patch = tensor.extract_image_patches(1, 1);
  VERIFY_IS_EQUAL(single_pixel_patch.dimension(0), 2);
  VERIFY_IS_EQUAL(single_pixel_patch.dimension(1), 1);
  VERIFY_IS_EQUAL(single_pixel_patch.dimension(2), 1);
  VERIFY_IS_EQUAL(single_pixel_patch.dimension(3), 3 * 5);
  VERIFY_IS_EQUAL(single_pixel_patch.dimension(4), 7);

  // Single pixel patch: RowMajor
  Tensor<float, 5, RowMajor> single_pixel_patch_row_major;
  single_pixel_patch_row_major = tensor_row_major.extract_image_patches(1, 1);
  VERIFY_IS_EQUAL(single_pixel_patch_row_major.dimension(0), 7);
  VERIFY_IS_EQUAL(single_pixel_patch_row_major.dimension(1), 3 * 5);
  VERIFY_IS_EQUAL(single_pixel_patch_row_major.dimension(2), 1);
  VERIFY_IS_EQUAL(single_pixel_patch_row_major.dimension(3), 1);
  VERIFY_IS_EQUAL(single_pixel_patch_row_major.dimension(4), 2);

  for (int i = 0; i < tensor.size(); ++i) {
    // ColMajor
    if (tensor.data()[i] != single_pixel_patch.data()[i]) {
      std::cout << "Mismatch detected at index " << i << " : " << tensor.data()[i] << " vs "
                << single_pixel_patch.data()[i] << std::endl;
    }
    VERIFY_IS_EQUAL(single_pixel_patch.data()[i], tensor.data()[i]);
    // RowMajor
    if (tensor_row_major.data()[i] != single_pixel_patch_row_major.data()[i]) {
      std::cout << "Mismatch detected at index " << i << " : " << tensor.data()[i] << " vs "
                << single_pixel_patch_row_major.data()[i] << std::endl;
    }
    VERIFY_IS_EQUAL(single_pixel_patch_row_major.data()[i], tensor_row_major.data()[i]);
    VERIFY_IS_EQUAL(tensor.data()[i], tensor_row_major.data()[i]);
    VERIFY_IS_EQUAL(single_pixel_patch.data()[i], single_pixel_patch_row_major.data()[i]);
  }

  // Entire image patch: ColMajor
  Tensor<float, 5> entire_image_patch;
  entire_image_patch = tensor.extract_image_patches(3, 5);
  VERIFY_IS_EQUAL(entire_image_patch.dimension(0), 2);
  VERIFY_IS_EQUAL(entire_image_patch.dimension(1), 3);
  VERIFY_IS_EQUAL(entire_image_patch.dimension(2), 5);
  VERIFY_IS_EQUAL(entire_image_patch.dimension(3), 3 * 5);
  VERIFY_IS_EQUAL(entire_image_patch.dimension(4), 7);

  // Entire image patch: RowMajor
  Tensor<float, 5, RowMajor> entire_image_patch_row_major;
  entire_image_patch_row_major = tensor_row_major.extract_image_patches(3, 5);
  VERIFY_IS_EQUAL(entire_image_patch_row_major.dimension(0), 7);
  VERIFY_IS_EQUAL(entire_image_patch_row_major.dimension(1), 3 * 5);
  VERIFY_IS_EQUAL(entire_image_patch_row_major.dimension(2), 5);
  VERIFY_IS_EQUAL(entire_image_patch_row_major.dimension(3), 3);
  VERIFY_IS_EQUAL(entire_image_patch_row_major.dimension(4), 2);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 5; ++j) {
      int patchId = i + 3 * j;
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 5; ++c) {
          for (int d = 0; d < 2; ++d) {
            for (int b = 0; b < 7; ++b) {
              float expected = 0.0f;
              float expected_row_major = 0.0f;
              if (r - 1 + i >= 0 && c - 2 + j >= 0 && r - 1 + i < 3 && c - 2 + j < 5) {
                expected = tensor(d, r - 1 + i, c - 2 + j, b);
                expected_row_major = tensor_row_major(b, c - 2 + j, r - 1 + i, d);
              }
              // ColMajor
              if (entire_image_patch(d, r, c, patchId, b) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(entire_image_patch(d, r, c, patchId, b), expected);
              // RowMajor
              if (entire_image_patch_row_major(b, patchId, c, r, d) != expected_row_major) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(entire_image_patch_row_major(b, patchId, c, r, d), expected_row_major);
              // Check that ColMajor and RowMajor agree.
              VERIFY_IS_EQUAL(expected, expected_row_major);
            }
          }
        }
      }
    }
  }

  // 2D patch: ColMajor
  Tensor<float, 5> twod_patch;
  twod_patch = tensor.extract_image_patches(2, 2);
  VERIFY_IS_EQUAL(twod_patch.dimension(0), 2);
  VERIFY_IS_EQUAL(twod_patch.dimension(1), 2);
  VERIFY_IS_EQUAL(twod_patch.dimension(2), 2);
  VERIFY_IS_EQUAL(twod_patch.dimension(3), 3 * 5);
  VERIFY_IS_EQUAL(twod_patch.dimension(4), 7);

  // 2D patch: RowMajor
  Tensor<float, 5, RowMajor> twod_patch_row_major;
  twod_patch_row_major = tensor_row_major.extract_image_patches(2, 2);
  VERIFY_IS_EQUAL(twod_patch_row_major.dimension(0), 7);
  VERIFY_IS_EQUAL(twod_patch_row_major.dimension(1), 3 * 5);
  VERIFY_IS_EQUAL(twod_patch_row_major.dimension(2), 2);
  VERIFY_IS_EQUAL(twod_patch_row_major.dimension(3), 2);
  VERIFY_IS_EQUAL(twod_patch_row_major.dimension(4), 2);

  // Based on the calculation described in TensorTraits.h, padding happens to be 0.
  int row_padding = 0;
  int col_padding = 0;
  int stride = 1;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 5; ++j) {
      int patchId = i + 3 * j;
      for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 2; ++c) {
          for (int d = 0; d < 2; ++d) {
            for (int b = 0; b < 7; ++b) {
              float expected = 0.0f;
              float expected_row_major = 0.0f;
              int row_offset = r * stride + i - row_padding;
              int col_offset = c * stride + j - col_padding;
              // ColMajor
              if (row_offset >= 0 && col_offset >= 0 && row_offset < tensor.dimension(1) &&
                  col_offset < tensor.dimension(2)) {
                expected = tensor(d, row_offset, col_offset, b);
              }
              if (twod_patch(d, r, c, patchId, b) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(twod_patch(d, r, c, patchId, b), expected);

              // RowMajor
              if (row_offset >= 0 && col_offset >= 0 && row_offset < tensor_row_major.dimension(2) &&
                  col_offset < tensor_row_major.dimension(1)) {
                expected_row_major = tensor_row_major(b, col_offset, row_offset, d);
              }
              if (twod_patch_row_major(b, patchId, c, r, d) != expected_row_major) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(twod_patch_row_major(b, patchId, c, r, d), expected_row_major);
              // Check that ColMajor and RowMajor agree.
              VERIFY_IS_EQUAL(expected, expected_row_major);
            }
          }
        }
      }
    }
  }
}

// Verifies VALID padding (no padding) with incrementing values.
void test_patch_padding_valid() {
  int input_depth = 3;
  int input_rows = 3;
  int input_cols = 3;
  int input_batches = 1;
  int ksize = 2;   // Corresponds to the Rows and Cols for tensor.extract_image_patches<>.
  int stride = 2;  // Only same stride is supported.
  Tensor<float, 4> tensor(input_depth, input_rows, input_cols, input_batches);
  // Initializes tensor with incrementing numbers.
  for (int i = 0; i < tensor.size(); ++i) {
    tensor.data()[i] = i + 1;
  }
  // ColMajor
  Tensor<float, 5> result = tensor.extract_image_patches(ksize, ksize, stride, stride, 1, 1, PADDING_VALID);

  VERIFY_IS_EQUAL(result.dimension(0), input_depth);    // depth
  VERIFY_IS_EQUAL(result.dimension(1), ksize);          // kernel rows
  VERIFY_IS_EQUAL(result.dimension(2), ksize);          // kernel cols
  VERIFY_IS_EQUAL(result.dimension(3), 1);              // number of patches
  VERIFY_IS_EQUAL(result.dimension(4), input_batches);  // number of batches

  // RowMajor
  Tensor<float, 4, RowMajor> tensor_row_major = tensor.swap_layout();
  VERIFY_IS_EQUAL(tensor.dimension(0), tensor_row_major.dimension(3));
  VERIFY_IS_EQUAL(tensor.dimension(1), tensor_row_major.dimension(2));
  VERIFY_IS_EQUAL(tensor.dimension(2), tensor_row_major.dimension(1));
  VERIFY_IS_EQUAL(tensor.dimension(3), tensor_row_major.dimension(0));

  Tensor<float, 5, RowMajor> result_row_major =
      tensor_row_major.extract_image_patches(ksize, ksize, stride, stride, 1, 1, PADDING_VALID);
  VERIFY_IS_EQUAL(result.dimension(0), result_row_major.dimension(4));
  VERIFY_IS_EQUAL(result.dimension(1), result_row_major.dimension(3));
  VERIFY_IS_EQUAL(result.dimension(2), result_row_major.dimension(2));
  VERIFY_IS_EQUAL(result.dimension(3), result_row_major.dimension(1));
  VERIFY_IS_EQUAL(result.dimension(4), result_row_major.dimension(0));

  // No padding is carried out.
  int row_padding = 0;
  int col_padding = 0;

  for (int i = 0; (i + stride + ksize - 1) < input_rows; i += stride) {    // input rows
    for (int j = 0; (j + stride + ksize - 1) < input_cols; j += stride) {  // input cols
      int patchId = i + input_rows * j;
      for (int r = 0; r < ksize; ++r) {                // patch rows
        for (int c = 0; c < ksize; ++c) {              // patch cols
          for (int d = 0; d < input_depth; ++d) {      // depth
            for (int b = 0; b < input_batches; ++b) {  // batch
              float expected = 0.0f;
              float expected_row_major = 0.0f;
              int row_offset = r + i - row_padding;
              int col_offset = c + j - col_padding;
              if (row_offset >= 0 && col_offset >= 0 && row_offset < input_rows && col_offset < input_cols) {
                expected = tensor(d, row_offset, col_offset, b);
                expected_row_major = tensor_row_major(b, col_offset, row_offset, d);
              }
              // ColMajor
              if (result(d, r, c, patchId, b) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(result(d, r, c, patchId, b), expected);
              // RowMajor
              if (result_row_major(b, patchId, c, r, d) != expected_row_major) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(result_row_major(b, patchId, c, r, d), expected_row_major);
              // Check that ColMajor and RowMajor agree.
              VERIFY_IS_EQUAL(expected, expected_row_major);
            }
          }
        }
      }
    }
  }
}

// Verifies VALID padding (no padding) with the same value.
void test_patch_padding_valid_same_value() {
  int input_depth = 1;
  int input_rows = 5;
  int input_cols = 5;
  int input_batches = 2;
  int ksize = 3;   // Corresponds to the Rows and Cols for tensor.extract_image_patches<>.
  int stride = 2;  // Only same stride is supported.
  // ColMajor
  Tensor<float, 4> tensor(input_depth, input_rows, input_cols, input_batches);
  tensor = tensor.constant(11.0f);
  Tensor<float, 5> result = tensor.extract_image_patches(ksize, ksize, stride, stride, 1, 1, PADDING_VALID);

  VERIFY_IS_EQUAL(result.dimension(0), input_depth);    // depth
  VERIFY_IS_EQUAL(result.dimension(1), ksize);          // kernel rows
  VERIFY_IS_EQUAL(result.dimension(2), ksize);          // kernel cols
  VERIFY_IS_EQUAL(result.dimension(3), 4);              // number of patches
  VERIFY_IS_EQUAL(result.dimension(4), input_batches);  // number of batches

  // RowMajor
  Tensor<float, 4, RowMajor> tensor_row_major = tensor.swap_layout();
  VERIFY_IS_EQUAL(tensor.dimension(0), tensor_row_major.dimension(3));
  VERIFY_IS_EQUAL(tensor.dimension(1), tensor_row_major.dimension(2));
  VERIFY_IS_EQUAL(tensor.dimension(2), tensor_row_major.dimension(1));
  VERIFY_IS_EQUAL(tensor.dimension(3), tensor_row_major.dimension(0));

  Tensor<float, 5, RowMajor> result_row_major =
      tensor_row_major.extract_image_patches(ksize, ksize, stride, stride, 1, 1, PADDING_VALID);
  VERIFY_IS_EQUAL(result.dimension(0), result_row_major.dimension(4));
  VERIFY_IS_EQUAL(result.dimension(1), result_row_major.dimension(3));
  VERIFY_IS_EQUAL(result.dimension(2), result_row_major.dimension(2));
  VERIFY_IS_EQUAL(result.dimension(3), result_row_major.dimension(1));
  VERIFY_IS_EQUAL(result.dimension(4), result_row_major.dimension(0));

  // No padding is carried out.
  int row_padding = 0;
  int col_padding = 0;

  for (int i = 0; (i + stride + ksize - 1) <= input_rows; i += stride) {    // input rows
    for (int j = 0; (j + stride + ksize - 1) <= input_cols; j += stride) {  // input cols
      int patchId = i + input_rows * j;
      for (int r = 0; r < ksize; ++r) {                // patch rows
        for (int c = 0; c < ksize; ++c) {              // patch cols
          for (int d = 0; d < input_depth; ++d) {      // depth
            for (int b = 0; b < input_batches; ++b) {  // batch
              float expected = 0.0f;
              float expected_row_major = 0.0f;
              int row_offset = r + i - row_padding;
              int col_offset = c + j - col_padding;
              if (row_offset >= 0 && col_offset >= 0 && row_offset < input_rows && col_offset < input_cols) {
                expected = tensor(d, row_offset, col_offset, b);
                expected_row_major = tensor_row_major(b, col_offset, row_offset, d);
              }
              // ColMajor
              if (result(d, r, c, patchId, b) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(result(d, r, c, patchId, b), expected);
              // RowMajor
              if (result_row_major(b, patchId, c, r, d) != expected_row_major) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(result_row_major(b, patchId, c, r, d), expected_row_major);
              // Check that ColMajor and RowMajor agree.
              VERIFY_IS_EQUAL(expected, expected_row_major);
            }
          }
        }
      }
    }
  }
}

// Verifies SAME padding.
void test_patch_padding_same() {
  int input_depth = 3;
  int input_rows = 4;
  int input_cols = 2;
  int input_batches = 1;
  int ksize = 2;   // Corresponds to the Rows and Cols for tensor.extract_image_patches<>.
  int stride = 2;  // Only same stride is supported.
  // ColMajor
  Tensor<float, 4> tensor(input_depth, input_rows, input_cols, input_batches);
  // Initializes tensor with incrementing numbers.
  for (int i = 0; i < tensor.size(); ++i) {
    tensor.data()[i] = i + 1;
  }
  Tensor<float, 5> result = tensor.extract_image_patches(ksize, ksize, stride, stride, PADDING_SAME);

  VERIFY_IS_EQUAL(result.dimension(0), input_depth);    // depth
  VERIFY_IS_EQUAL(result.dimension(1), ksize);          // kernel rows
  VERIFY_IS_EQUAL(result.dimension(2), ksize);          // kernel cols
  VERIFY_IS_EQUAL(result.dimension(3), 2);              // number of patches
  VERIFY_IS_EQUAL(result.dimension(4), input_batches);  // number of batches

  // RowMajor
  Tensor<float, 4, RowMajor> tensor_row_major = tensor.swap_layout();
  VERIFY_IS_EQUAL(tensor.dimension(0), tensor_row_major.dimension(3));
  VERIFY_IS_EQUAL(tensor.dimension(1), tensor_row_major.dimension(2));
  VERIFY_IS_EQUAL(tensor.dimension(2), tensor_row_major.dimension(1));
  VERIFY_IS_EQUAL(tensor.dimension(3), tensor_row_major.dimension(0));

  Tensor<float, 5, RowMajor> result_row_major =
      tensor_row_major.extract_image_patches(ksize, ksize, stride, stride, PADDING_SAME);
  VERIFY_IS_EQUAL(result.dimension(0), result_row_major.dimension(4));
  VERIFY_IS_EQUAL(result.dimension(1), result_row_major.dimension(3));
  VERIFY_IS_EQUAL(result.dimension(2), result_row_major.dimension(2));
  VERIFY_IS_EQUAL(result.dimension(3), result_row_major.dimension(1));
  VERIFY_IS_EQUAL(result.dimension(4), result_row_major.dimension(0));

  // Based on the calculation described in TensorTraits.h, padding happens to be
  // 0.
  int row_padding = 0;
  int col_padding = 0;

  for (int i = 0; (i + stride + ksize - 1) <= input_rows; i += stride) {    // input rows
    for (int j = 0; (j + stride + ksize - 1) <= input_cols; j += stride) {  // input cols
      int patchId = i + input_rows * j;
      for (int r = 0; r < ksize; ++r) {                // patch rows
        for (int c = 0; c < ksize; ++c) {              // patch cols
          for (int d = 0; d < input_depth; ++d) {      // depth
            for (int b = 0; b < input_batches; ++b) {  // batch
              float expected = 0.0f;
              float expected_row_major = 0.0f;
              int row_offset = r * stride + i - row_padding;
              int col_offset = c * stride + j - col_padding;
              if (row_offset >= 0 && col_offset >= 0 && row_offset < input_rows && col_offset < input_cols) {
                expected = tensor(d, row_offset, col_offset, b);
                expected_row_major = tensor_row_major(b, col_offset, row_offset, d);
              }
              // ColMajor
              if (result(d, r, c, patchId, b) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(result(d, r, c, patchId, b), expected);
              // RowMajor
              if (result_row_major(b, patchId, c, r, d) != expected_row_major) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(result_row_major(b, patchId, c, r, d), expected_row_major);
              // Check that ColMajor and RowMajor agree.
              VERIFY_IS_EQUAL(expected, expected_row_major);
            }
          }
        }
      }
    }
  }
}

// Verifies that SAME padding, when computed as negative values, will be clipped
// to zero.
void test_patch_padding_same_negative_padding_clip_to_zero() {
  int input_depth = 1;
  int input_rows = 15;
  int input_cols = 1;
  int input_batches = 1;
  int ksize = 1;  // Corresponds to the Rows and Cols for
                  // tensor.extract_image_patches<>.
  int row_stride = 5;
  int col_stride = 1;
  // ColMajor
  Tensor<float, 4> tensor(input_depth, input_rows, input_cols, input_batches);
  // Initializes tensor with incrementing numbers.
  for (int i = 0; i < tensor.size(); ++i) {
    tensor.data()[i] = i + 1;
  }
  Tensor<float, 5> result = tensor.extract_image_patches(ksize, ksize, row_stride, col_stride, 1, 1, PADDING_SAME);
  // row padding will be computed as -2 originally and then be clipped to 0.
  VERIFY_IS_EQUAL(result.coeff(0), 1.0f);
  VERIFY_IS_EQUAL(result.coeff(1), 6.0f);
  VERIFY_IS_EQUAL(result.coeff(2), 11.0f);

  VERIFY_IS_EQUAL(result.dimension(0), input_depth);    // depth
  VERIFY_IS_EQUAL(result.dimension(1), ksize);          // kernel rows
  VERIFY_IS_EQUAL(result.dimension(2), ksize);          // kernel cols
  VERIFY_IS_EQUAL(result.dimension(3), 3);              // number of patches
  VERIFY_IS_EQUAL(result.dimension(4), input_batches);  // number of batches

  // RowMajor
  Tensor<float, 4, RowMajor> tensor_row_major = tensor.swap_layout();
  VERIFY_IS_EQUAL(tensor.dimension(0), tensor_row_major.dimension(3));
  VERIFY_IS_EQUAL(tensor.dimension(1), tensor_row_major.dimension(2));
  VERIFY_IS_EQUAL(tensor.dimension(2), tensor_row_major.dimension(1));
  VERIFY_IS_EQUAL(tensor.dimension(3), tensor_row_major.dimension(0));

  Tensor<float, 5, RowMajor> result_row_major =
      tensor_row_major.extract_image_patches(ksize, ksize, row_stride, col_stride, 1, 1, PADDING_SAME);
  VERIFY_IS_EQUAL(result_row_major.coeff(0), 1.0f);
  VERIFY_IS_EQUAL(result_row_major.coeff(1), 6.0f);
  VERIFY_IS_EQUAL(result_row_major.coeff(2), 11.0f);

  VERIFY_IS_EQUAL(result.dimension(0), result_row_major.dimension(4));
  VERIFY_IS_EQUAL(result.dimension(1), result_row_major.dimension(3));
  VERIFY_IS_EQUAL(result.dimension(2), result_row_major.dimension(2));
  VERIFY_IS_EQUAL(result.dimension(3), result_row_major.dimension(1));
  VERIFY_IS_EQUAL(result.dimension(4), result_row_major.dimension(0));
}

void test_patch_no_extra_dim() {
  Tensor<float, 3> tensor(2, 3, 5);
  tensor.setRandom();
  Tensor<float, 3, RowMajor> tensor_row_major = tensor.swap_layout();
  VERIFY_IS_EQUAL(tensor.dimension(0), tensor_row_major.dimension(2));
  VERIFY_IS_EQUAL(tensor.dimension(1), tensor_row_major.dimension(1));
  VERIFY_IS_EQUAL(tensor.dimension(2), tensor_row_major.dimension(0));

  // Single pixel patch: ColMajor
  Tensor<float, 4> single_pixel_patch;
  single_pixel_patch = tensor.extract_image_patches(1, 1);
  VERIFY_IS_EQUAL(single_pixel_patch.dimension(0), 2);
  VERIFY_IS_EQUAL(single_pixel_patch.dimension(1), 1);
  VERIFY_IS_EQUAL(single_pixel_patch.dimension(2), 1);
  VERIFY_IS_EQUAL(single_pixel_patch.dimension(3), 3 * 5);

  // Single pixel patch: RowMajor
  Tensor<float, 4, RowMajor> single_pixel_patch_row_major;
  single_pixel_patch_row_major = tensor_row_major.extract_image_patches(1, 1);
  VERIFY_IS_EQUAL(single_pixel_patch_row_major.dimension(0), 3 * 5);
  VERIFY_IS_EQUAL(single_pixel_patch_row_major.dimension(1), 1);
  VERIFY_IS_EQUAL(single_pixel_patch_row_major.dimension(2), 1);
  VERIFY_IS_EQUAL(single_pixel_patch_row_major.dimension(3), 2);

  for (int i = 0; i < tensor.size(); ++i) {
    // ColMajor
    if (tensor.data()[i] != single_pixel_patch.data()[i]) {
      std::cout << "Mismatch detected at index " << i << " : " << tensor.data()[i] << " vs "
                << single_pixel_patch.data()[i] << std::endl;
    }
    VERIFY_IS_EQUAL(single_pixel_patch.data()[i], tensor.data()[i]);
    // RowMajor
    if (tensor_row_major.data()[i] != single_pixel_patch_row_major.data()[i]) {
      std::cout << "Mismatch detected at index " << i << " : " << tensor.data()[i] << " vs "
                << single_pixel_patch_row_major.data()[i] << std::endl;
    }
    VERIFY_IS_EQUAL(single_pixel_patch_row_major.data()[i], tensor_row_major.data()[i]);
    VERIFY_IS_EQUAL(tensor.data()[i], tensor_row_major.data()[i]);
    VERIFY_IS_EQUAL(single_pixel_patch.data()[i], single_pixel_patch_row_major.data()[i]);
  }

  // Entire image patch: ColMajor
  Tensor<float, 4> entire_image_patch;
  entire_image_patch = tensor.extract_image_patches(3, 5);
  VERIFY_IS_EQUAL(entire_image_patch.dimension(0), 2);
  VERIFY_IS_EQUAL(entire_image_patch.dimension(1), 3);
  VERIFY_IS_EQUAL(entire_image_patch.dimension(2), 5);
  VERIFY_IS_EQUAL(entire_image_patch.dimension(3), 3 * 5);

  // Entire image patch: RowMajor
  Tensor<float, 4, RowMajor> entire_image_patch_row_major;
  entire_image_patch_row_major = tensor_row_major.extract_image_patches(3, 5);
  VERIFY_IS_EQUAL(entire_image_patch_row_major.dimension(0), 3 * 5);
  VERIFY_IS_EQUAL(entire_image_patch_row_major.dimension(1), 5);
  VERIFY_IS_EQUAL(entire_image_patch_row_major.dimension(2), 3);
  VERIFY_IS_EQUAL(entire_image_patch_row_major.dimension(3), 2);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 5; ++j) {
      int patchId = i + 3 * j;
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 5; ++c) {
          for (int d = 0; d < 2; ++d) {
            float expected = 0.0f;
            float expected_row_major = 0.0f;
            if (r - 1 + i >= 0 && c - 2 + j >= 0 && r - 1 + i < 3 && c - 2 + j < 5) {
              expected = tensor(d, r - 1 + i, c - 2 + j);
              expected_row_major = tensor_row_major(c - 2 + j, r - 1 + i, d);
            }
            // ColMajor
            if (entire_image_patch(d, r, c, patchId) != expected) {
              std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c << " d=" << d
                        << std::endl;
            }
            VERIFY_IS_EQUAL(entire_image_patch(d, r, c, patchId), expected);
            // RowMajor
            if (entire_image_patch_row_major(patchId, c, r, d) != expected_row_major) {
              std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c << " d=" << d
                        << std::endl;
            }
            VERIFY_IS_EQUAL(entire_image_patch_row_major(patchId, c, r, d), expected_row_major);
            // Check that ColMajor and RowMajor agree.
            VERIFY_IS_EQUAL(expected, expected_row_major);
          }
        }
      }
    }
  }

  // 2D patch: ColMajor
  Tensor<float, 4> twod_patch;
  twod_patch = tensor.extract_image_patches(2, 2);
  VERIFY_IS_EQUAL(twod_patch.dimension(0), 2);
  VERIFY_IS_EQUAL(twod_patch.dimension(1), 2);
  VERIFY_IS_EQUAL(twod_patch.dimension(2), 2);
  VERIFY_IS_EQUAL(twod_patch.dimension(3), 3 * 5);

  // 2D patch: RowMajor
  Tensor<float, 4, RowMajor> twod_patch_row_major;
  twod_patch_row_major = tensor_row_major.extract_image_patches(2, 2);
  VERIFY_IS_EQUAL(twod_patch_row_major.dimension(0), 3 * 5);
  VERIFY_IS_EQUAL(twod_patch_row_major.dimension(1), 2);
  VERIFY_IS_EQUAL(twod_patch_row_major.dimension(2), 2);
  VERIFY_IS_EQUAL(twod_patch_row_major.dimension(3), 2);

  // Based on the calculation described in TensorTraits.h, padding happens to be 0.
  int row_padding = 0;
  int col_padding = 0;
  int stride = 1;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 5; ++j) {
      int patchId = i + 3 * j;
      for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 2; ++c) {
          for (int d = 0; d < 2; ++d) {
            float expected = 0.0f;
            float expected_row_major = 0.0f;
            int row_offset = r * stride + i - row_padding;
            int col_offset = c * stride + j - col_padding;
            // ColMajor
            if (row_offset >= 0 && col_offset >= 0 && row_offset < tensor.dimension(1) &&
                col_offset < tensor.dimension(2)) {
              expected = tensor(d, row_offset, col_offset);
            }
            if (twod_patch(d, r, c, patchId) != expected) {
              std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c << " d=" << d
                        << std::endl;
            }
            VERIFY_IS_EQUAL(twod_patch(d, r, c, patchId), expected);
            // RowMajor
            if (row_offset >= 0 && col_offset >= 0 && row_offset < tensor_row_major.dimension(1) &&
                col_offset < tensor_row_major.dimension(0)) {
              expected_row_major = tensor_row_major(col_offset, row_offset, d);
            }
            if (twod_patch_row_major(patchId, c, r, d) != expected_row_major) {
              std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c << " d=" << d
                        << std::endl;
            }
            VERIFY_IS_EQUAL(twod_patch_row_major(patchId, c, r, d), expected_row_major);
            // Check that ColMajor and RowMajor agree.
            VERIFY_IS_EQUAL(expected, expected_row_major);
          }
        }
      }
    }
  }
}

void test_imagenet_patches() {
  // Test the code on typical configurations used by the 'imagenet' benchmarks at
  // https://github.com/soumith/convnet-benchmarks
  // ColMajor
  Tensor<float, 4> l_in(3, 128, 128, 16);
  l_in.setRandom();
  Tensor<float, 5> l_out = l_in.extract_image_patches(11, 11);
  VERIFY_IS_EQUAL(l_out.dimension(0), 3);
  VERIFY_IS_EQUAL(l_out.dimension(1), 11);
  VERIFY_IS_EQUAL(l_out.dimension(2), 11);
  VERIFY_IS_EQUAL(l_out.dimension(3), 128 * 128);
  VERIFY_IS_EQUAL(l_out.dimension(4), 16);

  // RowMajor
  Tensor<float, 5, RowMajor> l_out_row_major = l_in.swap_layout().extract_image_patches(11, 11);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(0), 16);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(1), 128 * 128);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(2), 11);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(3), 11);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(4), 3);

  for (int b = 0; b < 16; ++b) {
    for (int i = 0; i < 128; ++i) {
      for (int j = 0; j < 128; ++j) {
        int patchId = i + 128 * j;
        for (int c = 0; c < 11; ++c) {
          for (int r = 0; r < 11; ++r) {
            for (int d = 0; d < 3; ++d) {
              float expected = 0.0f;
              if (r - 5 + i >= 0 && c - 5 + j >= 0 && r - 5 + i < 128 && c - 5 + j < 128) {
                expected = l_in(d, r - 5 + i, c - 5 + j, b);
              }
              // ColMajor
              if (l_out(d, r, c, patchId, b) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(l_out(d, r, c, patchId, b), expected);
              // RowMajor
              if (l_out_row_major(b, patchId, c, r, d) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(l_out_row_major(b, patchId, c, r, d), expected);
            }
          }
        }
      }
    }
  }

  // ColMajor
  l_in.resize(16, 64, 64, 32);
  l_in.setRandom();
  l_out = l_in.extract_image_patches(9, 9);
  VERIFY_IS_EQUAL(l_out.dimension(0), 16);
  VERIFY_IS_EQUAL(l_out.dimension(1), 9);
  VERIFY_IS_EQUAL(l_out.dimension(2), 9);
  VERIFY_IS_EQUAL(l_out.dimension(3), 64 * 64);
  VERIFY_IS_EQUAL(l_out.dimension(4), 32);

  // RowMajor
  l_out_row_major = l_in.swap_layout().extract_image_patches(9, 9);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(0), 32);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(1), 64 * 64);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(2), 9);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(3), 9);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(4), 16);

  for (int b = 0; b < 32; ++b) {
    for (int i = 0; i < 64; ++i) {
      for (int j = 0; j < 64; ++j) {
        int patchId = i + 64 * j;
        for (int c = 0; c < 9; ++c) {
          for (int r = 0; r < 9; ++r) {
            for (int d = 0; d < 16; ++d) {
              float expected = 0.0f;
              if (r - 4 + i >= 0 && c - 4 + j >= 0 && r - 4 + i < 64 && c - 4 + j < 64) {
                expected = l_in(d, r - 4 + i, c - 4 + j, b);
              }
              // ColMajor
              if (l_out(d, r, c, patchId, b) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(l_out(d, r, c, patchId, b), expected);
              // RowMajor
              if (l_out_row_major(b, patchId, c, r, d) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(l_out_row_major(b, patchId, c, r, d), expected);
            }
          }
        }
      }
    }
  }

  // ColMajor
  l_in.resize(32, 16, 16, 32);
  l_in.setRandom();
  l_out = l_in.extract_image_patches(7, 7);
  VERIFY_IS_EQUAL(l_out.dimension(0), 32);
  VERIFY_IS_EQUAL(l_out.dimension(1), 7);
  VERIFY_IS_EQUAL(l_out.dimension(2), 7);
  VERIFY_IS_EQUAL(l_out.dimension(3), 16 * 16);
  VERIFY_IS_EQUAL(l_out.dimension(4), 32);

  // RowMajor
  l_out_row_major = l_in.swap_layout().extract_image_patches(7, 7);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(0), 32);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(1), 16 * 16);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(2), 7);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(3), 7);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(4), 32);

  for (int b = 0; b < 32; ++b) {
    for (int i = 0; i < 16; ++i) {
      for (int j = 0; j < 16; ++j) {
        int patchId = i + 16 * j;
        for (int c = 0; c < 7; ++c) {
          for (int r = 0; r < 7; ++r) {
            for (int d = 0; d < 32; ++d) {
              float expected = 0.0f;
              if (r - 3 + i >= 0 && c - 3 + j >= 0 && r - 3 + i < 16 && c - 3 + j < 16) {
                expected = l_in(d, r - 3 + i, c - 3 + j, b);
              }
              // ColMajor
              if (l_out(d, r, c, patchId, b) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(l_out(d, r, c, patchId, b), expected);
              // RowMajor
              if (l_out_row_major(b, patchId, c, r, d) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(l_out_row_major(b, patchId, c, r, d), expected);
            }
          }
        }
      }
    }
  }

  // ColMajor
  l_in.resize(64, 13, 13, 32);
  l_in.setRandom();
  l_out = l_in.extract_image_patches(3, 3);
  VERIFY_IS_EQUAL(l_out.dimension(0), 64);
  VERIFY_IS_EQUAL(l_out.dimension(1), 3);
  VERIFY_IS_EQUAL(l_out.dimension(2), 3);
  VERIFY_IS_EQUAL(l_out.dimension(3), 13 * 13);
  VERIFY_IS_EQUAL(l_out.dimension(4), 32);

  // RowMajor
  l_out_row_major = l_in.swap_layout().extract_image_patches(3, 3);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(0), 32);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(1), 13 * 13);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(2), 3);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(3), 3);
  VERIFY_IS_EQUAL(l_out_row_major.dimension(4), 64);

  for (int b = 0; b < 32; ++b) {
    for (int i = 0; i < 13; ++i) {
      for (int j = 0; j < 13; ++j) {
        int patchId = i + 13 * j;
        for (int c = 0; c < 3; ++c) {
          for (int r = 0; r < 3; ++r) {
            for (int d = 0; d < 64; ++d) {
              float expected = 0.0f;
              if (r - 1 + i >= 0 && c - 1 + j >= 0 && r - 1 + i < 13 && c - 1 + j < 13) {
                expected = l_in(d, r - 1 + i, c - 1 + j, b);
              }
              // ColMajor
              if (l_out(d, r, c, patchId, b) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(l_out(d, r, c, patchId, b), expected);
              // RowMajor
              if (l_out_row_major(b, patchId, c, r, d) != expected) {
                std::cout << "Mismatch detected at index i=" << i << " j=" << j << " r=" << r << " c=" << c
                          << " d=" << d << " b=" << b << std::endl;
              }
              VERIFY_IS_EQUAL(l_out_row_major(b, patchId, c, r, d), expected);
            }
          }
        }
      }
    }
  }
}

// Tests inflate strides (row_inflate_strides, col_inflate_strides).
// Inflate inserts zeros between input elements before patch extraction.
void test_patch_inflate_strides() {
  // ColMajor: 2 channels, 3 rows, 4 cols, 1 batch
  const int depth = 2;
  const int rows = 3;
  const int cols = 4;
  const int batch = 1;
  Tensor<float, 4> tensor(depth, rows, cols, batch);
  for (int i = 0; i < tensor.size(); ++i) {
    tensor.data()[i] = static_cast<float>(i + 1);
  }

  const int row_inflate = 2;
  const int col_inflate = 3;
  const int patch_rows = 3;
  const int patch_cols = 4;

  // Effective input size after inflation:
  //   eff_rows = (3-1)*2 + 1 = 5
  //   eff_cols = (4-1)*3 + 1 = 10
  // With explicit padding=0, VALID-like extraction:
  //   outputRows = ceil((5 + 0 + 0 - 3 + 1) / 1) = 3
  //   outputCols = ceil((10 + 0 + 0 - 4 + 1) / 1) = 7
  Tensor<float, 5> result =
      tensor.extract_image_patches(patch_rows, patch_cols, 1, 1, 1, 1, row_inflate, col_inflate, 0, 0, 0, 0, 0.0f);
  const int outputRows = 3;
  const int outputCols = 7;
  VERIFY_IS_EQUAL(result.dimension(0), depth);
  VERIFY_IS_EQUAL(result.dimension(1), patch_rows);
  VERIFY_IS_EQUAL(result.dimension(2), patch_cols);
  VERIFY_IS_EQUAL(result.dimension(3), outputRows * outputCols);
  VERIFY_IS_EQUAL(result.dimension(4), batch);

  for (int b = 0; b < batch; ++b) {
    for (int oi = 0; oi < outputRows; ++oi) {
      for (int oj = 0; oj < outputCols; ++oj) {
        int patchId = oi + outputRows * oj;
        for (int pr = 0; pr < patch_rows; ++pr) {
          for (int pc = 0; pc < patch_cols; ++pc) {
            for (int d = 0; d < depth; ++d) {
              // Position in effective (inflated) input
              int effRow = oi + pr;
              int effCol = oj + pc;
              float expected = 0.0f;
              // Check if this maps to an actual input element
              if (effRow % row_inflate == 0 && effCol % col_inflate == 0) {
                int origRow = effRow / row_inflate;
                int origCol = effCol / col_inflate;
                if (origRow >= 0 && origRow < rows && origCol >= 0 && origCol < cols) {
                  expected = tensor(d, origRow, origCol, b);
                }
              }
              VERIFY_IS_EQUAL(result(d, pr, pc, patchId, b), expected);
            }
          }
        }
      }
    }
  }

  // RowMajor
  Tensor<float, 4, RowMajor> tensor_rm = tensor.swap_layout();
  Tensor<float, 5, RowMajor> result_rm =
      tensor_rm.extract_image_patches(patch_rows, patch_cols, 1, 1, 1, 1, row_inflate, col_inflate, 0, 0, 0, 0, 0.0f);
  VERIFY_IS_EQUAL(result_rm.dimension(4), depth);
  VERIFY_IS_EQUAL(result_rm.dimension(3), patch_rows);
  VERIFY_IS_EQUAL(result_rm.dimension(2), patch_cols);
  VERIFY_IS_EQUAL(result_rm.dimension(1), outputRows * outputCols);
  VERIFY_IS_EQUAL(result_rm.dimension(0), batch);

  for (int b = 0; b < batch; ++b) {
    for (int oi = 0; oi < outputRows; ++oi) {
      for (int oj = 0; oj < outputCols; ++oj) {
        int patchId = oi + outputRows * oj;
        for (int pr = 0; pr < patch_rows; ++pr) {
          for (int pc = 0; pc < patch_cols; ++pc) {
            for (int d = 0; d < depth; ++d) {
              int effRow = oi + pr;
              int effCol = oj + pc;
              float expected = 0.0f;
              if (effRow % row_inflate == 0 && effCol % col_inflate == 0) {
                int origRow = effRow / row_inflate;
                int origCol = effCol / col_inflate;
                if (origRow >= 0 && origRow < rows && origCol >= 0 && origCol < cols) {
                  expected = tensor_rm(b, origCol, origRow, d);
                }
              }
              VERIFY_IS_EQUAL(result_rm(b, patchId, pc, pr, d), expected);
            }
          }
        }
      }
    }
  }
}

// Tests dilation (in_row_strides, in_col_strides).
// Dilation samples every Nth element within each patch.
void test_patch_dilation() {
  const int depth = 3;
  const int rows = 5;
  const int cols = 5;
  const int batch = 1;
  Tensor<float, 4> tensor(depth, rows, cols, batch);
  for (int i = 0; i < tensor.size(); ++i) {
    tensor.data()[i] = static_cast<float>(i + 1);
  }

  const int patch_rows = 2;
  const int patch_cols = 2;
  const int in_row_strides = 2;  // dilation
  const int in_col_strides = 2;

  // Effective patch size: patch + (patch-1)*(dilation-1)
  //   eff_patch_rows = 2 + (2-1)*(2-1) = 3
  //   eff_patch_cols = 2 + (2-1)*(2-1) = 3
  // With PADDING_VALID:
  //   outputRows = ceil((5 - 3 + 1) / 1) = 3
  //   outputCols = ceil((5 - 3 + 1) / 1) = 3
  Tensor<float, 5> result =
      tensor.extract_image_patches(patch_rows, patch_cols, 1, 1, in_row_strides, in_col_strides, PADDING_VALID);
  const int outputRows = 3;
  const int outputCols = 3;
  VERIFY_IS_EQUAL(result.dimension(0), depth);
  VERIFY_IS_EQUAL(result.dimension(1), patch_rows);
  VERIFY_IS_EQUAL(result.dimension(2), patch_cols);
  VERIFY_IS_EQUAL(result.dimension(3), outputRows * outputCols);
  VERIFY_IS_EQUAL(result.dimension(4), batch);

  // row_padding and col_padding are 0 for VALID.
  for (int b = 0; b < batch; ++b) {
    for (int oi = 0; oi < outputRows; ++oi) {
      for (int oj = 0; oj < outputCols; ++oj) {
        int patchId = oi + outputRows * oj;
        for (int pr = 0; pr < patch_rows; ++pr) {
          for (int pc = 0; pc < patch_cols; ++pc) {
            for (int d = 0; d < depth; ++d) {
              // Within-patch dilation: sample at stride in_row_strides
              int inputRow = oi + pr * in_row_strides;
              int inputCol = oj + pc * in_col_strides;
              float expected = 0.0f;
              if (inputRow >= 0 && inputRow < rows && inputCol >= 0 && inputCol < cols) {
                expected = tensor(d, inputRow, inputCol, b);
              }
              VERIFY_IS_EQUAL(result(d, pr, pc, patchId, b), expected);
            }
          }
        }
      }
    }
  }

  // RowMajor
  Tensor<float, 4, RowMajor> tensor_rm = tensor.swap_layout();
  Tensor<float, 5, RowMajor> result_rm =
      tensor_rm.extract_image_patches(patch_rows, patch_cols, 1, 1, in_row_strides, in_col_strides, PADDING_VALID);
  VERIFY_IS_EQUAL(result_rm.dimension(4), depth);
  VERIFY_IS_EQUAL(result_rm.dimension(3), patch_rows);
  VERIFY_IS_EQUAL(result_rm.dimension(2), patch_cols);
  VERIFY_IS_EQUAL(result_rm.dimension(1), outputRows * outputCols);
  VERIFY_IS_EQUAL(result_rm.dimension(0), batch);

  for (int b = 0; b < batch; ++b) {
    for (int oi = 0; oi < outputRows; ++oi) {
      for (int oj = 0; oj < outputCols; ++oj) {
        int patchId = oi + outputRows * oj;
        for (int pr = 0; pr < patch_rows; ++pr) {
          for (int pc = 0; pc < patch_cols; ++pc) {
            for (int d = 0; d < depth; ++d) {
              int inputRow = oi + pr * in_row_strides;
              int inputCol = oj + pc * in_col_strides;
              float expected = 0.0f;
              if (inputRow >= 0 && inputRow < rows && inputCol >= 0 && inputCol < cols) {
                expected = tensor_rm(b, inputCol, inputRow, d);
              }
              VERIFY_IS_EQUAL(result_rm(b, patchId, pc, pr, d), expected);
            }
          }
        }
      }
    }
  }
}

// Tests explicit padding with asymmetric top/bottom/left/right values.
void test_patch_explicit_padding() {
  const int depth = 3;
  const int rows = 4;
  const int cols = 4;
  const int batch = 1;
  Tensor<float, 4> tensor(depth, rows, cols, batch);
  for (int i = 0; i < tensor.size(); ++i) {
    tensor.data()[i] = static_cast<float>(i + 1);
  }

  const int patch_rows = 3;
  const int patch_cols = 3;
  const int padding_top = 1;
  const int padding_bottom = 2;
  const int padding_left = 1;
  const int padding_right = 2;

  // outputRows = ceil((4 + 1 + 2 - 3 + 1) / 1) = 5
  // outputCols = ceil((4 + 1 + 2 - 3 + 1) / 1) = 5
  Tensor<float, 5> result = tensor.extract_image_patches(patch_rows, patch_cols, 1, 1, 1, 1, 1, 1, padding_top,
                                                         padding_bottom, padding_left, padding_right, 0.0f);
  const int outputRows = 5;
  const int outputCols = 5;
  VERIFY_IS_EQUAL(result.dimension(0), depth);
  VERIFY_IS_EQUAL(result.dimension(1), patch_rows);
  VERIFY_IS_EQUAL(result.dimension(2), patch_cols);
  VERIFY_IS_EQUAL(result.dimension(3), outputRows * outputCols);
  VERIFY_IS_EQUAL(result.dimension(4), batch);

  for (int b = 0; b < batch; ++b) {
    for (int oi = 0; oi < outputRows; ++oi) {
      for (int oj = 0; oj < outputCols; ++oj) {
        int patchId = oi + outputRows * oj;
        for (int pr = 0; pr < patch_rows; ++pr) {
          for (int pc = 0; pc < patch_cols; ++pc) {
            for (int d = 0; d < depth; ++d) {
              int inputRow = oi + pr - padding_top;
              int inputCol = oj + pc - padding_left;
              float expected = 0.0f;
              if (inputRow >= 0 && inputRow < rows && inputCol >= 0 && inputCol < cols) {
                expected = tensor(d, inputRow, inputCol, b);
              }
              VERIFY_IS_EQUAL(result(d, pr, pc, patchId, b), expected);
            }
          }
        }
      }
    }
  }

  // RowMajor
  Tensor<float, 4, RowMajor> tensor_rm = tensor.swap_layout();
  Tensor<float, 5, RowMajor> result_rm = tensor_rm.extract_image_patches(
      patch_rows, patch_cols, 1, 1, 1, 1, 1, 1, padding_top, padding_bottom, padding_left, padding_right, 0.0f);
  VERIFY_IS_EQUAL(result_rm.dimension(4), depth);
  VERIFY_IS_EQUAL(result_rm.dimension(3), patch_rows);
  VERIFY_IS_EQUAL(result_rm.dimension(2), patch_cols);
  VERIFY_IS_EQUAL(result_rm.dimension(1), outputRows * outputCols);
  VERIFY_IS_EQUAL(result_rm.dimension(0), batch);

  for (int b = 0; b < batch; ++b) {
    for (int oi = 0; oi < outputRows; ++oi) {
      for (int oj = 0; oj < outputCols; ++oj) {
        int patchId = oi + outputRows * oj;
        for (int pr = 0; pr < patch_rows; ++pr) {
          for (int pc = 0; pc < patch_cols; ++pc) {
            for (int d = 0; d < depth; ++d) {
              int inputRow = oi + pr - padding_top;
              int inputCol = oj + pc - padding_left;
              float expected = 0.0f;
              if (inputRow >= 0 && inputRow < rows && inputCol >= 0 && inputCol < cols) {
                expected = tensor_rm(b, inputCol, inputRow, d);
              }
              VERIFY_IS_EQUAL(result_rm(b, patchId, pc, pr, d), expected);
            }
          }
        }
      }
    }
  }
}

// Tests rectangular input with non-square patches and different row/col strides.
void test_patch_asymmetric() {
  const int depth = 3;
  const int rows = 3;
  const int cols = 7;
  const int batch = 1;
  Tensor<float, 4> tensor(depth, rows, cols, batch);
  for (int i = 0; i < tensor.size(); ++i) {
    tensor.data()[i] = static_cast<float>(i + 1);
  }

  const int patch_rows = 2;
  const int patch_cols = 3;
  const int row_stride = 1;
  const int col_stride = 2;

  // PADDING_VALID:
  //   outputRows = ceil((3 - 2 + 1) / 1) = 2
  //   outputCols = ceil((7 - 3 + 1) / 2) = 3
  Tensor<float, 5> result =
      tensor.extract_image_patches(patch_rows, patch_cols, row_stride, col_stride, 1, 1, PADDING_VALID);
  const int outputRows = 2;
  const int outputCols = 3;
  VERIFY_IS_EQUAL(result.dimension(0), depth);
  VERIFY_IS_EQUAL(result.dimension(1), patch_rows);
  VERIFY_IS_EQUAL(result.dimension(2), patch_cols);
  VERIFY_IS_EQUAL(result.dimension(3), outputRows * outputCols);
  VERIFY_IS_EQUAL(result.dimension(4), batch);

  for (int b = 0; b < batch; ++b) {
    for (int oi = 0; oi < outputRows; ++oi) {
      for (int oj = 0; oj < outputCols; ++oj) {
        int patchId = oi + outputRows * oj;
        for (int pr = 0; pr < patch_rows; ++pr) {
          for (int pc = 0; pc < patch_cols; ++pc) {
            for (int d = 0; d < depth; ++d) {
              int inputRow = oi * row_stride + pr;
              int inputCol = oj * col_stride + pc;
              float expected = 0.0f;
              if (inputRow >= 0 && inputRow < rows && inputCol >= 0 && inputCol < cols) {
                expected = tensor(d, inputRow, inputCol, b);
              }
              VERIFY_IS_EQUAL(result(d, pr, pc, patchId, b), expected);
            }
          }
        }
      }
    }
  }

  // RowMajor
  Tensor<float, 4, RowMajor> tensor_rm = tensor.swap_layout();
  Tensor<float, 5, RowMajor> result_rm =
      tensor_rm.extract_image_patches(patch_rows, patch_cols, row_stride, col_stride, 1, 1, PADDING_VALID);
  VERIFY_IS_EQUAL(result_rm.dimension(4), depth);
  VERIFY_IS_EQUAL(result_rm.dimension(3), patch_rows);
  VERIFY_IS_EQUAL(result_rm.dimension(2), patch_cols);
  VERIFY_IS_EQUAL(result_rm.dimension(1), outputRows * outputCols);
  VERIFY_IS_EQUAL(result_rm.dimension(0), batch);

  for (int b = 0; b < batch; ++b) {
    for (int oi = 0; oi < outputRows; ++oi) {
      for (int oj = 0; oj < outputCols; ++oj) {
        int patchId = oi + outputRows * oj;
        for (int pr = 0; pr < patch_rows; ++pr) {
          for (int pc = 0; pc < patch_cols; ++pc) {
            for (int d = 0; d < depth; ++d) {
              int inputRow = oi * row_stride + pr;
              int inputCol = oj * col_stride + pc;
              float expected = 0.0f;
              if (inputRow >= 0 && inputRow < rows && inputCol >= 0 && inputCol < cols) {
                expected = tensor_rm(b, inputCol, inputRow, d);
              }
              VERIFY_IS_EQUAL(result_rm(b, patchId, pc, pr, d), expected);
            }
          }
        }
      }
    }
  }
}

// Exercises packet loads that span multiple rows/columns within the patch.
void test_patch_contiguous_packet_span() {
  const int depth = 1;
  const int rows = 5;
  const int cols = 5;
  const int batch = 1;
  Tensor<float, 4> tensor(depth, rows, cols, batch);
  for (int i = 0; i < tensor.size(); ++i) {
    tensor.data()[i] = static_cast<float>(i + 1);
  }

  const int patch_rows = 3;
  const int patch_cols = 3;
  Tensor<float, 5> result = tensor.extract_image_patches(patch_rows, patch_cols, 1, 1, 1, 1, PADDING_VALID);
  const int outputRows = 3;
  const int outputCols = 3;

  for (int b = 0; b < batch; ++b) {
    for (int oi = 0; oi < outputRows; ++oi) {
      for (int oj = 0; oj < outputCols; ++oj) {
        const int patchId = oi + outputRows * oj;
        for (int pr = 0; pr < patch_rows; ++pr) {
          for (int pc = 0; pc < patch_cols; ++pc) {
            const int inputRow = oi + pr;
            const int inputCol = oj + pc;
            VERIFY_IS_EQUAL(result(0, pr, pc, patchId, b), tensor(0, inputRow, inputCol, b));
          }
        }
      }
    }
  }

  Tensor<float, 4, RowMajor> tensor_rm = tensor.swap_layout();
  Tensor<float, 5, RowMajor> result_rm =
      tensor_rm.extract_image_patches(patch_rows, patch_cols, 1, 1, 1, 1, PADDING_VALID);
  for (int b = 0; b < batch; ++b) {
    for (int oi = 0; oi < outputRows; ++oi) {
      for (int oj = 0; oj < outputCols; ++oj) {
        const int patchId = oi + outputRows * oj;
        for (int pr = 0; pr < patch_rows; ++pr) {
          for (int pc = 0; pc < patch_cols; ++pc) {
            const int inputRow = oi + pr;
            const int inputCol = oj + pc;
            VERIFY_IS_EQUAL(result_rm(b, patchId, pc, pr, 0), tensor_rm(b, inputCol, inputRow, 0));
          }
        }
      }
    }
  }
}

EIGEN_DECLARE_TEST(tensor_image_patch) {
  CALL_SUBTEST_1(test_simple_patch());
  CALL_SUBTEST_2(test_patch_no_extra_dim());
  CALL_SUBTEST_3(test_patch_padding_valid());
  CALL_SUBTEST_4(test_patch_padding_valid_same_value());
  CALL_SUBTEST_5(test_patch_padding_same());
  CALL_SUBTEST_6(test_imagenet_patches());
  CALL_SUBTEST_7(test_patch_padding_same_negative_padding_clip_to_zero());
  CALL_SUBTEST_8(test_patch_inflate_strides());
  CALL_SUBTEST_9(test_patch_dilation());
  CALL_SUBTEST_10(test_patch_explicit_padding());
  CALL_SUBTEST_11(test_patch_asymmetric());
  CALL_SUBTEST_12(test_patch_contiguous_packet_span());
}
