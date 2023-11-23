#pragma once

#ifndef UTILS_SHA3_HH
#define UTILS_SHA3_HH

#include "Utils.hh"

#include <string>
#include <iostream>
#include <cstdint>

namespace Utils {

  /// SHA-3 winning hash algorithm Keccak
  ///
  /// @author: Christopher Bentivenga
  /// @author: Frederick Christie
  /// @author: Michael Kitson
  ///
  /// Porting in C++ by Enrico Bertolazzi
  ///
  class SHA3 {

    using string       = std::string;
    using ostream_type = std::basic_ostream<char>;
    using istream_type = std::basic_istream<char>;
    using uint8_t      = std::uint8_t;
    using uint64_t     = std::uint64_t;

    // Round state
    uint8_t  *m_buffer_location; // used for writing and to know when to flush the buffer
    uint64_t m_state[5][5];
    uint64_t m_message_buffer_64[1600/8];  // rate bits wide, defined during construction

    int m_digest_size; // bytes

    // Digest-length specific Values
    int m_sponge_capacity;
    int m_sponge_rate;

    void m_reset();
    void m_perform_rounds( int rounds );
    void m_absorb_buffer();

    // Debugging
    void m_print_message_buffer( ostream_type & stream ) const;
    void m_print_sponge( ostream_type & stream ) const;

  public:

    using keccakLane_t = uint64_t;

    SHA3( int digest_size );
    ~SHA3() { }

    /// Adds an entire string to the message
    ///
    /// @param  str  The string of bytes to add
    void hash_string( char const * str );

    /// Adds an entire hexidecimal string to the message
    ///
    /// @param  str  The hex string of bytes to add
    void hash_hex_string( char const * str );

    /// Returns a representation of the digest as a hexidecimal string
    ///
    /// @return The hex string, ownership of which is given to the caller
    string digest_in_hex();

    // Overridden functions from HashFunction
    int digest_size() const { return m_digest_size; }

    void hash( const int b );

    void digest( uint8_t * d );
  };
}

#endif
