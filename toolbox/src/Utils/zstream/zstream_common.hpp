/*
 zstream-cpp Library License:
 --------------------------

 The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

 This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.

 Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:

 1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.

 2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.

 3. This notice may not be removed or altered from any source distribution

 Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003
         Gero Mueller, post@geromueller.de, 2015
*/

#ifndef ZIP_STREAM_COMMON_HPP
#define ZIP_STREAM_COMMON_HPP

#include <stdint.h>

namespace zstream {

/// Compression strategy, see zlib doc.
enum EStrategy {
	StrategyFiltered = 1, StrategyHuffmanOnly = 2, DefaultStrategy = 0
};

namespace detail {
/// default gzip buffer size,
/// change this to suite your needs
const size_t default_buffer_size = 4096;

const int gz_magic[2] = { 0x1f, 0x8b }; /* gzip magic header */

/* gzip flag byte */
const int gz_ascii_flag = 0x01; /* bit 0 set: file probably ascii text */
const int gz_head_crc = 0x02; /* bit 1 set: header CRC present */
const int gz_extra_field = 0x04; /* bit 2 set: extra field present */
const int gz_orig_name = 0x08; /* bit 3 set: original file name present */
const int gz_comment = 0x10; /* bit 4 set: file comment present */
const int gz_reserved = 0xE0; /* bits 5..7: reserved */

} // detail
} // zstream

#endif

