/*
 zstream-cpp Library License:
 --------------------------

 The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

 This software is provided 'as-is', without any express or
 implied warranty. In no event will the authors be held
 liable for any damages arising from the use of this software.

 Permission is granted to anyone to use this software
 for any purpose, including commercial applications,
 and to alter it and redistribute it freely, subject
 to the following restrictions:

 1. The origin of this software must not be misrepresented;
    you must not claim that you wrote the original software.
    If you use this software in a product, an acknowledgment
    in the product documentation would be appreciated but is not required.

 2. Altered source versions must be plainly marked as such,
    and must not be misrepresented as being the original software.

 3. This notice may not be removed or altered from any source distribution

 Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003
         Gero Mueller, post@geromueller.de, 2015
*/

#ifndef OUTPUT_ZIP_STREAM_IMPL_HPP
#define OUTPUT_ZIP_STREAM_IMPL_HPP

#include "ozstream.hpp"

#ifndef OS_CODE
#  define OS_CODE  0x03  /* assume Unix */
#endif

namespace zstream {

template <
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::basic_zip_streambuf(
  ostream_reference ostream_,
  size_t            level_,
  EStrategy         strategy_,
  size_t            window_size_,
  size_t            memory_level_,
  size_t            buffer_size_
) 
: m_ostream(ostream_)
, m_output_buffer(buffer_size_, 0)
, m_buffer(buffer_size_, 0)
, m_crc(0)
{
  m_zip_stream.zalloc = (alloc_func) 0;
  m_zip_stream.zfree = (free_func) 0;

  m_zip_stream.next_in = NULL;
  m_zip_stream.avail_in = 0;
  m_zip_stream.avail_out = 0;
  m_zip_stream.next_out = NULL;

  m_err = deflateInit2(
	&m_zip_stream, std::min(9, static_cast<int>(level_)),
    Z_DEFLATED,
    -static_cast<int>(window_size_), // <-- changed
    std::min(9, static_cast<int>(memory_level_)),
	static_cast<int>(strategy_)
  );

  m_zip_stream.avail_out = static_cast<uInt>(m_output_buffer.size());
  m_zip_stream.next_out  = &(m_output_buffer[0]);

  char_type *p = &(m_buffer[0]);
  this->setp(p, p + buffer_size_);
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::~basic_zip_streambuf() {
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
int
basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::sync() {
  int c = overflow(EOF);
  if ( c == EOF ) return -1;
  else            return 0;
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
typename basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::int_type
basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::overflow(
  typename basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::int_type c
) {
  // buffer full, zip it..
  int w = static_cast<int>(this->pptr() - this->pbase());
  if (w > 0) {
    if (zip_to_stream(this->pbase(), w)) {
      this->setp(this->pbase(), this->epptr());
    } else {
      return EOF;
    }
  }

  // write c if not EOF
  if (c != EOF) {
    *this->pptr() = c;
    this->pbump(1);
    return c;
  } else {
    return 0;
  }
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
std::streamsize
basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::zip_write( int flag ) {
  std::streamsize total_written_byte_size = 0;

  m_err = ::deflate(&m_zip_stream, flag);

  if ( m_err == Z_OK || m_err == Z_STREAM_END ) {
	std::streamsize written_byte_size =
	  static_cast<std::streamsize>(m_output_buffer.size()) - m_zip_stream.avail_out;
    if (written_byte_size > 0) {
      total_written_byte_size += written_byte_size;
      // ouput buffer is full, dumping to ostream
      m_ostream.write(
		    (const char_type*) &(m_output_buffer[0]),
        static_cast<std::streamsize>(written_byte_size/sizeof(char_type))
      );

      // checking if some bytes were not written.
      size_t remainder = written_byte_size % sizeof(char_type);
	    if ( remainder != 0 ) {
	      // copy to the beginning of the stream
        memcpy(
          &(m_output_buffer[0]),
          &(m_output_buffer[written_byte_size - remainder]),
          remainder
        );
      }

      m_zip_stream.avail_out = static_cast<uInt>(m_output_buffer.size() - remainder);
      m_zip_stream.next_out = &m_output_buffer[remainder];
    }
  }
  return total_written_byte_size;
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
bool
basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::zip_to_stream(
  typename basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::char_type* buffer_,
  std::streamsize buffer_size_
) {
  m_zip_stream.next_in  = (byte_buffer_type) buffer_;
  m_zip_stream.avail_in = static_cast<uInt>(buffer_size_ * sizeof(char_type));

  // updating crc
  m_crc = ::crc32(m_crc, m_zip_stream.next_in, m_zip_stream.avail_in);

  do {
    zip_write(0);
  } while (m_zip_stream.avail_in != 0 && m_err == Z_OK);

  return m_err == Z_OK;
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
std::streamsize
basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::zfinish() {
  std::streamsize total_written_byte_size = 0;

  this->sync();

  // updating crc
  m_crc = crc32(m_crc, m_zip_stream.next_in, m_zip_stream.avail_in);

  do {
    total_written_byte_size += zip_write(Z_FINISH);
  } while (m_err == Z_OK);

  m_ostream.flush();
  m_err = deflateEnd(&m_zip_stream);

  return total_written_byte_size;
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
void
basic_gzip_ostream<Elem, Tr, ElemA, ByteT, ByteAT>::put_long(
  typename basic_gzip_ostream<Elem, Tr, ElemA, ByteT, ByteAT>::ostream_reference out_,
  unsigned int x_
) {
  static const int size_ul = sizeof(unsigned int);
  static const int size_c = sizeof(char_type);
  static const int n_end = size_ul / size_c;
  out_.write(reinterpret_cast<char_type const*>(&x_), n_end);
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
void
basic_gzip_ostream<Elem, Tr, ElemA, ByteT, ByteAT>::add_header() {
  char_type zero = 0;
  this->rdbuf()->get_ostream().put(
    static_cast<char_type>(detail::gz_magic[0])).put(
    static_cast<char_type>(detail::gz_magic[1])).put(
    static_cast<char_type>(Z_DEFLATED)
  )
  .put(zero) //flags
  .put(zero).put(zero).put(zero).put(zero) // time
  .put(zero) //xflags
  .put(static_cast<char_type>(OS_CODE));
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
void
basic_gzip_ostream<Elem, Tr, ElemA, ByteT, ByteAT>::add_footer() {
  put_long(this->rdbuf()->get_ostream(), this->rdbuf()->get_crc());
  put_long(this->rdbuf()->get_ostream(), this->rdbuf()->get_in_size());
}

} // zstream

#endif
