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

#ifndef OUTPUT_ZIP_STREAM_HPP
#define OUTPUT_ZIP_STREAM_HPP

#include "zstream_common.hpp"

namespace zstream {

//!
//! A stream decorator that takes raw input and zips it to a ostream.
//! 
//! The class wraps up the inflate method of the zlib library
//!
template<
  typename Elem,
  typename Tr = std::char_traits<Elem>,
  typename ElemA = std::allocator<Elem>,
  typename ByteT = unsigned char,
  typename ByteAT = std::allocator<ByteT>
>
class basic_zip_streambuf: public std::basic_streambuf<Elem, Tr> {
public:
  typedef std::basic_ostream<Elem, Tr>& ostream_reference;
  typedef ElemA char_allocator_type;
  typedef ByteT byte_type;
  typedef ByteAT byte_allocator_type;
  typedef byte_type* byte_buffer_type;
  typedef typename std::basic_streambuf<Elem, Tr>::char_type char_type;
  typedef typename std::basic_streambuf<Elem, Tr>::int_type int_type;
  typedef std::vector<byte_type, byte_allocator_type> byte_vector_type;
  typedef std::vector<char_type, char_allocator_type> char_vector_type;

  //! Construct a zip stream.
  //!
  //! More info on the following parameters can be
  //! found in the zlib documentation.
  //! 
  basic_zip_streambuf(
    ostream_reference ostream_,
    size_t            level_,
    EStrategy         strategy_,
    size_t            window_size_,
    size_t            memory_level_,
    size_t            buffer_size_
  );

  ~basic_zip_streambuf() override;

  int sync() override;
  int_type overflow( int_type c ) override;

  //!
  //! Flushes the zip buffer and output buffer.
  //! 
  //! This method should be called at the end of the compression. Calling flush multiple times, will lower the
  //! compression ratio.
  //!
  std::streamsize zfinish();

  //! returns a reference to the output stream
  ostream_reference get_ostream() const { return m_ostream; }

  //! returns the latest zlib error status
  int get_zerr() const { return m_err; }

  //! returns the crc of the input data compressed so far.
  long get_crc() const { return m_crc; }

  //! returns the size (bytes) of the input data compressed so far.
  long get_in_size() const { return m_zip_stream.total_in; }

  //! returns the size (bytes) of the compressed data so far.
  long get_out_size() const { return m_zip_stream.total_out; }

private:
  bool zip_to_stream(char_type*, std::streamsize);
  std::streamsize zip_write(int flag);
  
  ostream_reference m_ostream;
  z_stream          m_zip_stream;
  int               m_err;
  byte_vector_type  m_output_buffer;
  char_vector_type  m_buffer;
  long              m_crc;
};

//!
//! Base class for zip ostreams.
//! 
//! Contains a basic_zip_streambuf.
//!
template<
  typename Elem,
  typename Tr = std::char_traits<Elem>,
  typename ElemA = std::allocator<Elem>,
  typename ByteT = unsigned char,
  typename ByteAT = std::allocator<ByteT>
>
class basic_zip_ostreambase: virtual public std::basic_ios<Elem, Tr> {
public:
  typedef std::basic_ostream<Elem, Tr>& ostream_reference;
  typedef basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT> zip_streambuf_type;
  
  //!
  //! Construct a zip stream.
  //!
  //! More info on the following parameters can be
  //! found in the zlib documentation.
  //! 
  basic_zip_ostreambase(
    ostream_reference ostream_,
    size_t            level_,
    EStrategy         strategy_,
    size_t            window_size_,
    size_t            memory_level_,
    size_t            buffer_size_
  ) : m_buf(ostream_, level_, strategy_, window_size_, memory_level_, buffer_size_) {
  	this->init(&m_buf);
  }

  //! Returns the underlying zip ostream object
  zip_streambuf_type* rdbuf() { return &m_buf; }

  //! Returns the zlib error state
  int get_zerr() const { return m_buf.get_zerr(); }

  //! Returns the uncompressed data crc
  long get_crc() const { return m_buf.get_crc(); }

  //! Returns the compressed data size
  long get_out_size() const { return m_buf.get_out_size(); }

  //! Returns the uncompressed data size
  long get_in_size() const { return m_buf.get_in_size(); }

private:
  zip_streambuf_type m_buf;
};

//!
//! A zipper ostream.
//!
//! This class is a ostream decorator that behaves 'almost'
//! like any other ostream.
//!
//! At construction, it takes any ostream that shall
//! be used to output of the compressed data.
//!
//! When finished, you need to call the special method
//! close or call the destructor
//! to flush all the intermidiate streams.
//!
//! Example:
//!
//! \code
//! // creating the target zip string, could be a fstream
//! ostringstream ostringstream_;
//! // creating the zip layer
//! zip_ostream zipper(ostringstream_);
//!
//!
//! // writing data
//! zipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
//! // zip ostream needs special flushing...
//! zipper.close();
//! \endcode
//!
template<
  typename Elem,
  typename Tr = std::char_traits<Elem>,
  typename ElemA = std::allocator<Elem>,
  typename ByteT = unsigned char,
  typename ByteAT = std::allocator<ByteT>
>
class basic_gzip_ostream
  : public basic_zip_ostreambase<Elem, Tr, ElemA, ByteT, ByteAT>
  , public std::basic_ostream<Elem, Tr> {
public:
  typedef typename std::basic_ostream<Elem, Tr>::char_type char_type;
  typedef std::basic_ostream<Elem, Tr>& ostream_reference;
  typedef basic_zip_ostreambase<Elem, Tr, ElemA, ByteT, ByteAT> zip_ostreambase_type;
  typedef std::basic_ostream<Elem, Tr> ostream_type;

  //!
  //! Constructs a zipper ostream decorator
  //!
  //! \param ostream_ ostream where the compressed output is written
  //! \param is_gzip_ true if gzip header and footer have to be added
  //! \param level_ level of compression 0, bad and fast, 9, good and slower,
  //! \param strategy_ compression strategy
  //! \param window_size_ see zlib doc
  //! \param memory_level_ see zlib doc
  //! \param buffer_size_ the buffer size used to zip data
  //!
  //! When is_gzip_ is true, a gzip header and footer is automatically added.
  //!
  basic_gzip_ostream(
    ostream_reference ostream_,
    size_t            level_        = Z_DEFAULT_COMPRESSION,
    EStrategy         strategy_     = DefaultStrategy, size_t window_size_ = 15,
    size_t            memory_level_ = 8,
	size_t            buffer_size_  = detail::default_buffer_size
  )
  : zip_ostreambase_type(ostream_, level_, strategy_, window_size_, memory_level_, buffer_size_)
  , ostream_type(this->rdbuf())
  ,	m_closed(false) {
    add_header();
  }

  ~basic_gzip_ostream() override { close(); }

  void
  close() {
    if ( m_closed ) return;
    this->flush();
    this->rdbuf()->zfinish();
    add_footer();
    m_closed = true;
  }

private:
  static void put_long(ostream_reference out_, unsigned int x_);
  void add_header();
  void add_footer();
  bool m_closed;
};

template<
  typename Elem,
  typename Tr     = std::char_traits<Elem>,
  typename ElemA  = std::allocator<Elem>,
  typename ByteT  = unsigned char,
  typename ByteAT = std::allocator<ByteT>
>
class basic_zip_ostream
  : public basic_zip_ostreambase<Elem, Tr, ElemA, ByteT, ByteAT>
  , public std::basic_ostream<Elem, Tr> {
public:
  typedef typename std::basic_ostream<Elem, Tr>::char_type char_type;
  typedef std::basic_ostream<Elem, Tr>& ostream_reference;
  typedef basic_zip_ostreambase<Elem, Tr, ElemA, ByteT, ByteAT> zip_ostreambase_type;
  typedef std::basic_ostream<Elem, Tr> ostream_type;
  
  //!
  //! Constructs a zipper ostream decorator
  //!
  //! \param ostream_ ostream where the compressed output is written
  //! \param is_gzip_ true if gzip header and footer have to be added
  //! \param level_ level of compression 0, bad and fast, 9, good and slower,
  //! \param strategy_ compression strategy
  //! \param window_size_ see zlib doc
  //! \param memory_level_ see zlib doc
  //! \param buffer_size_ the buffer size used to zip data
  //!
  //! When is_gzip_ is true, a gzip header and footer is automatically added.
  //!
  basic_zip_ostream(
    ostream_reference ostream_,
    size_t    level_        = Z_DEFAULT_COMPRESSION,
	EStrategy strategy_     = DefaultStrategy,
	size_t    window_size_  = 15,
	size_t    memory_level_ = 8,
	size_t    buffer_size_  = detail::default_buffer_size
  )
  : zip_ostreambase_type(ostream_, level_, strategy_, window_size_, memory_level_, buffer_size_)
  , ostream_type(this->rdbuf())
  , m_closed(false)
  { }

  ~basic_zip_ostream() { close(); }

  void
  close() {
    if ( m_closed ) return;
    this->flush();
    this->rdbuf()->zfinish();
    m_closed = true;
  }

private:
  bool m_closed;
};

//! A typedef for basic_zip_ostream<char>
typedef basic_gzip_ostream<char> ogzstream;
//! A typedef for basic_zip_ostream<wchar_t>
typedef basic_gzip_ostream<wchar_t> wgzostream;
//! A typedef for basic_zip_ostream<char>
typedef basic_zip_ostream<char> ozstream;
//! A typedef for basic_zip_ostream<wchar_t>
typedef basic_zip_ostream<wchar_t> wzostream;

} // zstream

#include "ozstream_impl.hpp"

#endif
