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

#ifndef INPUT_ZIP_STREAM_HPP
#define INPUT_ZIP_STREAM_HPP

#include "zstream_common.hpp"

namespace zstream {

//!
//! A stream decorator that takes compressed input and unzips it to a istream.
//! 
//! The class wraps up the deflate method of the zlib library 1.1.4 http://www.gzip.org/zlib/
//!
template<
  typename Elem,
  typename Tr = std::char_traits<Elem>,
  typename ElemA = std::allocator<Elem>,
  typename ByteT = unsigned char,
  typename ByteAT = std::allocator<ByteT>
>
class basic_unzip_streambuf: public std::basic_streambuf<Elem, Tr> {
public:
  typedef std::basic_istream<Elem, Tr>& istream_reference;
  typedef ElemA char_allocator_type;
  typedef ByteT byte_type;
  typedef ByteAT byte_allocator_type;
  typedef byte_type* byte_buffer_type;
  typedef typename std::basic_streambuf<Elem, Tr>::char_type char_type;
  typedef typename std::basic_streambuf<Elem, Tr>::int_type int_type;
  typedef std::vector<byte_type, byte_allocator_type> byte_vector_type;
  typedef std::vector<char_type, char_allocator_type> char_vector_type;

  //! Construct a unzip stream.
  //! More info on the following parameters can be found 
  //! in the zlib documentation.
  //! 
  basic_unzip_streambuf(
    istream_reference istream_,
    size_t            window_size_,
    size_t            read_buffer_size_,
    size_t            input_buffer_size_
  );

  ~basic_unzip_streambuf() override;

  int_type underflow() override;

  //! returns the compressed input istream
  istream_reference get_istream() { return m_istream; }

  //! returns the zlib stream structure
  z_stream& get_zip_stream() { return m_zip_stream; }

  //! returns the latest zlib error state
  int get_zerr() const { return m_err; }

  //! returns the crc of the uncompressed data so far 
  unsigned get_crc() const { return m_crc; }

  //! returns the number of uncompressed bytes
  unsigned get_out_size() const { return m_zip_stream.total_out; }

  //! returns the number of read compressed bytes
  unsigned long get_in_size() const { return m_zip_stream.total_in; }

private:

  void put_back_from_zip_stream();
  std::streamsize unzip_from_stream(char_type*, std::streamsize);

  size_t fill_input_buffer();

  istream_reference m_istream;
  z_stream          m_zip_stream;
  int               m_err;
  byte_vector_type  m_input_buffer;
  char_vector_type  m_buffer;
  long              m_crc;
};

//!
//! Base class for unzip istreams.
//! 
//! Contains a basic_unzip_streambuf.
//!
template<
  typename Elem,
  typename Tr = std::char_traits<Elem>,
  typename ElemA = std::allocator<Elem>,
  typename ByteT = unsigned char,
  typename ByteAT = std::allocator<ByteT>
>
class basic_zip_istreambase: virtual public std::basic_ios<Elem, Tr> {
public:
  typedef std::basic_istream<Elem, Tr>& istream_reference;
  typedef basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT> unzip_streambuf_type;

  basic_zip_istreambase(
    istream_reference ostream_,
    size_t            window_size_,
    size_t            read_buffer_size_,
    size_t            input_buffer_size_
  ) : m_buf(ostream_, window_size_, read_buffer_size_, input_buffer_size_) {
    this->init(&m_buf);
  }

  //! returns the underlying unzip istream object
  unzip_streambuf_type* rdbuf() { return &m_buf; }

  //! returns the zlib error state
  int get_zerr() const { return m_buf.get_zerr(); }

  //! returns the uncompressed data crc
  long get_crc() const { return m_buf.get_crc(); }

  //! returns the uncompressed data size
  long get_out_size() const { return m_buf.get_out_size(); }

  //! returns the compressed data size
  unsigned long get_in_size() const { return m_buf.get_in_size(); }

private:
  unzip_streambuf_type m_buf;
};

//!
//! A zipper istream.
//! 
//! This class is a istream decorator that behaves 'almost'
//! like any other ostream.
//! 
//! At construction, it takes any istream that shall 
//! be used to input of the compressed data.
//! 
//! Simlpe example:
//!
//! \code
//! // create a stream on zip string
//! istringstream istringstream_( ostringstream_.str());
//! // create unzipper istream
//! zip_istream unzipper( istringstream_);
//! 
//! // read and unzip
//! unzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;
//! \endcode
//!
template<
  typename Elem,
  typename Tr = std::char_traits<Elem>,
  typename ElemA = std::allocator<Elem>,
  typename ByteT = unsigned char,
  typename ByteAT = std::allocator<ByteT>
>
class basic_gzip_istream:
  public basic_zip_istreambase<Elem, Tr, ElemA, ByteT, ByteAT>,
  public std::basic_istream<Elem, Tr> {
public:
  typedef typename std::basic_istream<Elem, Tr>::char_type char_type;
  typedef std::basic_istream<Elem, Tr>& istream_reference;
  typedef basic_zip_istreambase<Elem, Tr, ElemA, ByteT, ByteAT> zip_istreambase_type;
  typedef std::basic_istream<Elem, Tr> istream_type;
  typedef unsigned char byte_type;

  //!
  //! Construct a unzipper stream
  //! 
  //! \param istream_ input buffer
  //! \param window_size_ 
  //! \param read_buffer_size_ 
  //! \param input_buffer_size_ 
  //! 
  basic_gzip_istream(
    istream_reference istream_,
    size_t window_size_ = 15,
    size_t read_buffer_size_ = detail::default_buffer_size,
    size_t input_buffer_size_ = detail::default_buffer_size
  )
  : zip_istreambase_type(istream_, window_size_, read_buffer_size_,	input_buffer_size_)
  , istream_type(this->rdbuf())
  , m_gzip_crc(0)
  , m_gzip_data_size(0)
  {
    if ( this->rdbuf()->get_zerr() == Z_OK ) check_header();
  }

  //! reads the gzip header
  void read_footer();
  //!
  //! Return crc check result
  //! 
  //! When you have finished reading the compressed data,
  //! call read_footer to read the uncompressed data crc.
  //! This method compares it to the crc of the uncompressed data.
  //! 
  //! \return true if crc check is succesful
  //!
  bool
  check_crc() const {
    return this->get_crc() == m_gzip_crc;
  }

  //! return data size check
  bool
  check_data_size() const {
    return this->get_out_size() == m_gzip_data_size;
  }

  //! return the crc value in the file
  unsigned get_gzip_crc() const { return m_gzip_crc; }

  //! return the data size in the file 
  unsigned get_gzip_data_size() const { return m_gzip_data_size; }

protected:
  static void read_long(istream_reference in_, unsigned int& x_);
  int check_header();
  unsigned m_gzip_crc;
  unsigned m_gzip_data_size;
};

template<
  typename Elem,
  typename Tr = std::char_traits<Elem>,
  typename ElemA = std::allocator<Elem>,
  typename ByteT = unsigned char,
  typename ByteAT = std::allocator<ByteT>
>
class basic_zip_istream
  : public basic_zip_istreambase<Elem, Tr, ElemA, ByteT,ByteAT>
  , public std::basic_istream<Elem, Tr> {
public:
  typedef typename std::basic_istream<Elem, Tr>::char_type char_type;
  typedef std::basic_istream<Elem, Tr>& istream_reference;
  typedef basic_zip_istreambase<Elem, Tr, ElemA, ByteT, ByteAT> zip_istreambase_type;
  typedef std::basic_istream<Elem, Tr> istream_type;
  typedef unsigned char byte_type;

  //!
  //! Construct a unzipper stream
  //! 
  //! \param istream_ input buffer
  //! \param window_size_
  //! \param read_buffer_size_
  //! \param input_buffer_size_
  //! 
  basic_zip_istream(
	istream_reference istream_,
    size_t window_size_ = 15,
    size_t read_buffer_size_ = detail::default_buffer_size,
    size_t input_buffer_size_ = detail::default_buffer_size
  )
  : zip_istreambase_type(istream_, window_size_, read_buffer_size_, input_buffer_size_)
  , istream_type(this->rdbuf())
  {}
};

//! A typedef for basic_zip_istream<char>
typedef basic_gzip_istream<char> igzstream;
//! A typedef for basic_zip_istream<wchart>
typedef basic_gzip_istream<wchar_t> wigzstream;

//! A typedef for basic_zip_istream<char>
typedef basic_zip_istream<char> izstream;
//! A typedef for basic_zip_istream<wchart>
typedef basic_zip_istream<wchar_t> wizstream;

} // zstream

#include "izstream_impl.hpp"

#endif

