/*
 zstream-cpp Library License:
 --------------------------

 The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

 This software is provided 'as-is', without any express 
 or implied warranty. In no event will the authors be held
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

#ifndef INPUT_ZIP_STREAM_IMPL_HPP
#define INPUT_ZIP_STREAM_IMPL_HPP

#include "izstream.hpp"

namespace zstream {

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::basic_unzip_streambuf(
  istream_reference istream_,
  size_t            window_size_,
  size_t            read_buffer_size_,
  size_t            input_buffer_size_
)
: m_istream(istream_)
, m_input_buffer(input_buffer_size_)
, m_buffer(read_buffer_size_)
, m_crc(0)
{
  // setting zalloc, zfree and opaque
  m_zip_stream.zalloc = (alloc_func) 0;
  m_zip_stream.zfree = (free_func) 0;
  
  m_zip_stream.next_in = NULL;
  m_zip_stream.avail_in = 0;
  m_zip_stream.avail_out = 0;
  m_zip_stream.next_out = NULL;

  m_err = inflateInit2(&m_zip_stream, -static_cast<int>(window_size_));

  this->setg(&(m_buffer[0]) + 4,     // beginning of putback area
             &(m_buffer[0]) + 4,     // read position
             &(m_buffer[0]) + 4);    // end position
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
size_t
basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::fill_input_buffer() {
  m_zip_stream.next_in = &(m_input_buffer[0]);
  m_istream.read(
    (char_type*) (&(m_input_buffer[0])),
    static_cast<std::streamsize>(m_input_buffer.size()/sizeof(char_type))
  );
  return m_zip_stream.avail_in = m_istream.gcount() * sizeof(char_type);
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::~basic_unzip_streambuf() {
  inflateEnd(&m_zip_stream);
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
typename basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::int_type
basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::underflow() {
  if (this->gptr() && (this->gptr() < this->egptr()))
    return *reinterpret_cast<unsigned char *>(this->gptr());

  int n_putback = static_cast<int>(this->gptr() - this->eback());
  if (n_putback > 4) n_putback = 4;

  memcpy(
    &(m_buffer[0]) + (4 - n_putback),
    this->gptr() - n_putback,
		n_putback * sizeof(char_type)
  );

  int num = unzip_from_stream(
	&(m_buffer[0]) + 4,
	static_cast<std::streamsize>((m_buffer.size() - 4) * sizeof(char_type))
  );
  if (num <= 0) // ERROR or EOF
    return EOF;

  // reset buffer pointers
  this->setg(
    &(m_buffer[0]) + (4 - n_putback),   // beginning of putback area
    &(m_buffer[0]) + 4,                 // read position
    &(m_buffer[0]) + 4 + num            // end of buffer
  );

  // return next character
  return *reinterpret_cast<unsigned char *>(this->gptr());
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
std::streamsize
basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::unzip_from_stream(
  char_type*      buffer_,
  std::streamsize buffer_size_
) {
  m_zip_stream.next_out  = (byte_buffer_type) buffer_;
  m_zip_stream.avail_out = static_cast<uInt>(buffer_size_ * sizeof(char_type));
  size_t count = m_zip_stream.avail_in;

  do {
    if (m_zip_stream.avail_in == 0) count = fill_input_buffer();
    if (m_zip_stream.avail_in) m_err = inflate(&m_zip_stream, Z_SYNC_FLUSH);
  } while (m_err == Z_OK && m_zip_stream.avail_out != 0 && count != 0);

  // updating crc
  m_crc = crc32(
  	m_crc,
  	(byte_buffer_type) buffer_,
  	buffer_size_ - m_zip_stream.avail_out / sizeof(char_type)
  );
  std::streamsize n_read = buffer_size_ - m_zip_stream.avail_out / sizeof(char_type);

  // check if it is the end
  if (m_err == Z_STREAM_END) put_back_from_zip_stream();
  return n_read;
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
void
basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::put_back_from_zip_stream() {
  if (m_zip_stream.avail_in == 0) return;

  m_istream.clear(std::ios::goodbit);
  m_istream.seekg(-static_cast<int>(m_zip_stream.avail_in),std::ios_base::cur);
  m_zip_stream.avail_in = 0;
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
int
basic_gzip_istream<Elem, Tr, ElemA, ByteT, ByteAT>::check_header() {
  int  method; /* method byte */
  int  flags;  /* flags byte */
  uInt len;
  int c;
  int err = 0;
  z_stream& zip_stream = this->rdbuf()->get_zip_stream();

  /* Check the gzip magic header */
  for (len = 0; len < 2; len++) {
    c = (int) this->rdbuf()->get_istream().get();
    if (c != detail::gz_magic[len]) {
      if (len != 0)
        this->rdbuf()->get_istream().unget();
      if (c != EOF) {
        this->rdbuf()->get_istream().unget();
      }
      err = zip_stream.avail_in != 0 ? Z_OK : Z_STREAM_END;
      return err;
	  }
  }

  method = (int) this->rdbuf()->get_istream().get();
  flags  = (int) this->rdbuf()->get_istream().get();
  if (method != Z_DEFLATED || (flags & detail::gz_reserved) != 0) {
    err = Z_DATA_ERROR;
    return err;
  }

  /* Discard time, xflags and OS code: */
  for (len = 0; len < 6; len++)
    this->rdbuf()->get_istream().get();

  if ((flags & detail::gz_extra_field) != 0) {
    /* skip the extra field */
    len = (uInt) this->rdbuf()->get_istream().get();
    len += ((uInt) this->rdbuf()->get_istream().get()) << 8;
    /* len is garbage if EOF but the loop below will quit anyway */
    while (len-- != 0 && this->rdbuf()->get_istream().get() != EOF)
      ;
  }
  if ((flags & detail::gz_orig_name) != 0) {
    /* skip the original file name */
    while ((c = this->rdbuf()->get_istream().get()) != 0 && c != EOF)
      ;
  }
  if ((flags & detail::gz_comment) != 0) {
    /* skip the .gz file comment */
    while ((c = this->rdbuf()->get_istream().get()) != 0 && c != EOF)
      ;
  }
  if ((flags & detail::gz_head_crc) != 0) { /* skip the header crc */
    for (len = 0; len < 2; len++)
      this->rdbuf()->get_istream().get();
  }
  err = this->rdbuf()->get_istream().eof() ? Z_DATA_ERROR : Z_OK;
  return err;
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
void
basic_gzip_istream<Elem, Tr, ElemA, ByteT, ByteAT>::read_footer() {
  read_long(this->rdbuf()->get_istream(), m_gzip_crc);
  read_long(this->rdbuf()->get_istream(), m_gzip_data_size);
}

template<
  typename Elem,
  typename Tr,
  typename ElemA,
  typename ByteT,
  typename ByteAT
>
void
basic_gzip_istream<Elem, Tr, ElemA, ByteT, ByteAT>::read_long(
  istream_reference in_,
  unsigned int& x_
) {
  static const int size_ul = sizeof(unsigned int);
  static const int size_c = sizeof(char_type);
  static const int n_end = size_ul / size_c;
  in_.read(reinterpret_cast<char*>(&x_), n_end);
}

} // zstream

#endif
