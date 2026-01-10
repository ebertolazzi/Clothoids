/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once

#ifndef UTILS_SHA3_HH
#define UTILS_SHA3_HH

#include "Utils.hh"

// Circular rotate left
#define ROT_L( X, Y ) ( ( ( X ) << ( Y ) ) | ( ( X ) >> ( 64 - ( Y ) ) ) )
#define ROUNDS 24

//! For converting binary output to hexadecimal for printing
static constexpr const char * hexLookup = "0123456789abcdef";

static constexpr uint64_t roundConstants[] = {
  0x0000000000000001, 0x0000000000008082, 0x800000000000808A, 0x8000000080008000, 0x000000000000808B,
  0x0000000080000001, 0x8000000080008081, 0x8000000000008009, 0x000000000000008A, 0x0000000000000088,
  0x0000000080008009, 0x000000008000000A, 0x000000008000808B, 0x800000000000008B, 0x8000000000008089,
  0x8000000000008003, 0x8000000000008002, 0x8000000000000080, 0x000000000000800A, 0x800000008000000A,
  0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008
};

namespace Utils
{

  //!
  //! \brief SHA-3 winning hash algorithm Keccak
  //!
  //! This class implements the SHA-3 hashing algorithm, based on the Keccak
  //! design. It provides methods to compute the hash of strings and returns
  //! the hash digest in hexadecimal format.
  //!
  //! \author Christopher Bentivenga
  //! \author Frederick Christie
  //! \author Michael Kitson
  //!
  //! Ported to C++ by Enrico Bertolazzi
  //!
  class SHA3
  {
    using uint8_t      = std::uint8_t;
    using uint64_t     = std::uint64_t;
    using ostream_type = std::ostream;
    using string_view  = const std::string &;

    // Round state
    uint8_t * m_buffer_location{ nullptr };  // used for writing and to know when to flush the buffer
    uint64_t  m_state[5][5];
    uint64_t  m_message_buffer_64[1600 / 8];  // rate bits wide, defined during
                                              // construction

    int m_digest_size{ 0 };  // bytes

    // Digest-length specific Values
    int m_sponge_capacity{ 0 };
    int m_sponge_rate{ 0 };

    // Track if we've already processed the final block
    bool m_finalized{ false };

    void m_reset()
    {
      auto const s{ static_cast<uint64_t *>( m_state[0] ) };
      std::fill_n( s, 25, 0 );
      m_buffer_location = reinterpret_cast<uint8_t *>( m_message_buffer_64 );
      m_finalized       = false;
    }

    void m_absorb_buffer()
    {
      uint64_t const * x{ m_message_buffer_64 };
      for ( int i{ 0 }; i * 64 < m_sponge_rate; ++i ) m_state[i / 5][i % 5] ^= x[i];  // TODO: unroll
      m_perform_rounds( ROUNDS );
    }

    // CHANGE: Function changed to inline
    void m_perform_rounds( int const rounds )
    {
      uint64_t b[5][5];
      uint64_t c[5];
      uint64_t d[5];

      for ( int i{ 0 }; i < rounds; i++ )
      {
        // CHANGE: For loops change to pre-determined steps, reduces branching

        // Theta step
        c[0] = m_state[0][0] ^ m_state[1][0] ^ m_state[2][0] ^ m_state[3][0] ^ m_state[4][0];
        c[1] = m_state[0][1] ^ m_state[1][1] ^ m_state[2][1] ^ m_state[3][1] ^ m_state[4][1];
        c[2] = m_state[0][2] ^ m_state[1][2] ^ m_state[2][2] ^ m_state[3][2] ^ m_state[4][2];
        c[3] = m_state[0][3] ^ m_state[1][3] ^ m_state[2][3] ^ m_state[3][3] ^ m_state[4][3];
        c[4] = m_state[0][4] ^ m_state[1][4] ^ m_state[2][4] ^ m_state[3][4] ^ m_state[4][4];

        d[0] = c[4] ^ ROT_L( c[1], 1 );
        d[1] = c[0] ^ ROT_L( c[2], 1 );
        d[2] = c[1] ^ ROT_L( c[3], 1 );
        d[3] = c[2] ^ ROT_L( c[4], 1 );
        d[4] = c[3] ^ ROT_L( c[0], 1 );

        m_state[0][0] ^= d[0];
        m_state[0][1] ^= d[1];
        m_state[0][2] ^= d[2];
        m_state[0][3] ^= d[3];
        m_state[0][4] ^= d[4];
        m_state[1][0] ^= d[0];
        m_state[1][1] ^= d[1];
        m_state[1][2] ^= d[2];
        m_state[1][3] ^= d[3];
        m_state[1][4] ^= d[4];
        m_state[2][0] ^= d[0];
        m_state[2][1] ^= d[1];
        m_state[2][2] ^= d[2];
        m_state[2][3] ^= d[3];
        m_state[2][4] ^= d[4];
        m_state[3][0] ^= d[0];
        m_state[3][1] ^= d[1];
        m_state[3][2] ^= d[2];
        m_state[3][3] ^= d[3];
        m_state[3][4] ^= d[4];
        m_state[4][0] ^= d[0];
        m_state[4][1] ^= d[1];
        m_state[4][2] ^= d[2];
        m_state[4][3] ^= d[3];
        m_state[4][4] ^= d[4];

        // Rho and Pi steps
        b[0][0] = m_state[0][0];  // rotate left by 0 bits
        b[1][3] = ROT_L( m_state[1][0], 36 );
        b[2][1] = ROT_L( m_state[2][0], 3 );
        b[3][4] = ROT_L( m_state[3][0], 41 );
        b[4][2] = ROT_L( m_state[4][0], 18 );

        b[0][2] = ROT_L( m_state[0][1], 1 );
        b[1][0] = ROT_L( m_state[1][1], 44 );
        b[2][3] = ROT_L( m_state[2][1], 10 );
        b[3][1] = ROT_L( m_state[3][1], 45 );
        b[4][4] = ROT_L( m_state[4][1], 2 );

        b[0][4] = ROT_L( m_state[0][2], 62 );
        b[1][2] = ROT_L( m_state[1][2], 6 );
        b[2][0] = ROT_L( m_state[2][2], 43 );
        b[3][3] = ROT_L( m_state[3][2], 15 );
        b[4][1] = ROT_L( m_state[4][2], 61 );

        b[0][1] = ROT_L( m_state[0][3], 28 );
        b[1][4] = ROT_L( m_state[1][3], 55 );
        b[2][2] = ROT_L( m_state[2][3], 25 );
        b[3][0] = ROT_L( m_state[3][3], 21 );
        b[4][3] = ROT_L( m_state[4][3], 56 );

        b[0][3] = ROT_L( m_state[0][4], 27 );
        b[1][1] = ROT_L( m_state[1][4], 20 );
        b[2][4] = ROT_L( m_state[2][4], 39 );
        b[3][2] = ROT_L( m_state[3][4], 8 );
        b[4][0] = ROT_L( m_state[4][4], 14 );

        // Chi step
        m_state[0][0] = b[0][0] ^ ( ( ~b[1][0] ) & b[2][0] );
        m_state[1][0] = b[0][1] ^ ( ( ~b[1][1] ) & b[2][1] );
        m_state[2][0] = b[0][2] ^ ( ( ~b[1][2] ) & b[2][2] );
        m_state[3][0] = b[0][3] ^ ( ( ~b[1][3] ) & b[2][3] );
        m_state[4][0] = b[0][4] ^ ( ( ~b[1][4] ) & b[2][4] );

        m_state[0][1] = b[1][0] ^ ( ( ~b[2][0] ) & b[3][0] );
        m_state[1][1] = b[1][1] ^ ( ( ~b[2][1] ) & b[3][1] );
        m_state[2][1] = b[1][2] ^ ( ( ~b[2][2] ) & b[3][2] );
        m_state[3][1] = b[1][3] ^ ( ( ~b[2][3] ) & b[3][3] );
        m_state[4][1] = b[1][4] ^ ( ( ~b[2][4] ) & b[3][4] );

        m_state[0][2] = b[2][0] ^ ( ( ~b[3][0] ) & b[4][0] );
        m_state[1][2] = b[2][1] ^ ( ( ~b[3][1] ) & b[4][1] );
        m_state[2][2] = b[2][2] ^ ( ( ~b[3][2] ) & b[4][2] );
        m_state[3][2] = b[2][3] ^ ( ( ~b[3][3] ) & b[4][3] );
        m_state[4][2] = b[2][4] ^ ( ( ~b[3][4] ) & b[4][4] );

        m_state[0][3] = b[3][0] ^ ( ( ~b[4][0] ) & b[0][0] );
        m_state[1][3] = b[3][1] ^ ( ( ~b[4][1] ) & b[0][1] );
        m_state[2][3] = b[3][2] ^ ( ( ~b[4][2] ) & b[0][2] );
        m_state[3][3] = b[3][3] ^ ( ( ~b[4][3] ) & b[0][3] );
        m_state[4][3] = b[3][4] ^ ( ( ~b[4][4] ) & b[0][4] );

        m_state[0][4] = b[4][0] ^ ( ( ~b[0][0] ) & b[1][0] );
        m_state[1][4] = b[4][1] ^ ( ( ~b[0][1] ) & b[1][1] );
        m_state[2][4] = b[4][2] ^ ( ( ~b[0][2] ) & b[1][2] );
        m_state[3][4] = b[4][3] ^ ( ( ~b[0][3] ) & b[1][3] );
        m_state[4][4] = b[4][4] ^ ( ( ~b[0][4] ) & b[1][4] );

        // Iota step
        m_state[0][0] ^= roundConstants[i];
      }
    }

    // Helper function for SHA3 padding
    void m_apply_sha3_padding()
    {
      // SHA3 padding: append 0x06 then pad with 0x80
      // This is equivalent to: append bits 01 (for SHA3) then apply Keccak padding 10*1

      uint8_t * buffer_start = reinterpret_cast<uint8_t *>( m_message_buffer_64 );
      uint8_t * buffer_end   = buffer_start + ( m_sponge_rate >> 3 );

      // Append 0x06 (binary: 00000110)
      // This represents: append 2 bits '01' (for SHA3 domain separation)
      // followed by first '1' of pad10*1
      *m_buffer_location = 0x06;
      m_buffer_location++;

      // Zero out the rest of the buffer
      if ( m_buffer_location < buffer_end ) { memset( m_buffer_location, 0, buffer_end - m_buffer_location ); }

      // Set the last byte to 0x80
      buffer_end[-1] = 0x80;

      m_absorb_buffer();
    }

    // Debugging
    void m_print_message_buffer( ostream_type & stream ) const
    {
      auto const m_messageBuffer{ reinterpret_cast<uint8_t const *>( m_message_buffer_64 ) };
      stream << "mb = [ ";
      for ( int i{ 0 }; i < m_sponge_rate / 8; ++i ) stream << static_cast<int>( m_messageBuffer[i] ) << ' ';
      stream << "]\n";
    }

    void m_print_sponge( ostream_type & stream ) const
    {
      stream << "s = [ " << std::hex;
      for ( int i{ 0 }; i < 5; ++i )
        for ( int j{ 0 }; j < 5; ++j ) stream << m_state[i][j] << ' ';
      stream << std::dec << "]\n";
    }

  public:
    //!
    //! \brief Constructor for the SHA3 class
    //!
    //! Initializes the SHA3 hashing algorithm with the specified digest size.
    //!
    //! \param digest_size The size of the desired hash output in bytes (e.g.,
    //! 224, 256, 384, or 512).
    //!
    explicit SHA3( int const digest_size ) : m_digest_size( digest_size )
    {
      // zero the state
      // CHANGE: Now uses bit shifting instead of multiplication
      m_sponge_capacity = m_digest_size << 4;  // x 16
      m_sponge_rate     = 1600 - m_sponge_capacity;
      m_reset();
    }

    //! \brief Destructor for the SHA3 class.
    ~SHA3() = default;

    //!
    //! \brief Hash a specified number of bytes.
    //!
    //! This method processes the provided byte input to update the hash state.
    //!
    //! \param b The number of bytes to hash.
    //!
    void hash( int const b )
    {
      if ( m_finalized ) { m_reset(); }

      m_buffer_location[0] = static_cast<uint8_t>( b );
      m_buffer_location++;
      if ( m_buffer_location == reinterpret_cast<uint8_t *>( m_message_buffer_64 ) + ( m_sponge_rate >> 3 ) )
      {
        m_buffer_location = reinterpret_cast<uint8_t *>( m_message_buffer_64 );
        m_absorb_buffer();
      }
    }

    //!
    //! \brief Adds an entire string to the message.
    //!
    //! This method appends the given string to the internal message buffer
    //! for hashing. It processes the string as a sequence of bytes.
    //!
    //! \param str The null-terminated string of bytes to add to the hash.
    //!
    void hash_string( string_view str )
    {
      if ( m_finalized ) { m_reset(); }

      int byte = 0;
      while ( str[byte] != '\0' )
      {
        hash( static_cast<int>( static_cast<uint8_t>( str[byte] ) ) );
        ++byte;
      }
    }

    //!
    //! \brief Adds an entire hexadecimal string to the message.
    //!
    //! This method appends the given hexadecimal string to the internal
    //! message buffer. Each pair of hex digits is converted to a byte.
    //!
    //! \param str The null-terminated hexadecimal string of bytes to add to the
    //! hash.
    //!
    void hash_hex_string( string_view str )
    {
      if ( m_finalized ) { m_reset(); }

      int byte{ 0 };
      while ( str[byte] != '\0' )
      {
        int f{ str[byte] };
        int s{ str[byte + 1] };
        if ( f >= 97 )
          f -= 87;  // lowercase
        else if ( f >= 65 )
          f -= 55;  // uppercase
        else
          f -= 48;  // numeric

        if ( s >= 97 )
          s -= 87;  // lowercase
        else if ( s >= 65 )
          s -= 55;  // uppercase
        else
          s -= 48;  // numeric

        hash( ( f << 4 ) | s );
        byte += 2;
      }
    }

    //!
    //! \brief Retrieve the digest and store it in the provided byte array.
    //!
    //! This method finalizes the hashing process and stores the resulting
    //! hash digest in the provided byte array.
    //!
    //! \param d Pointer to the byte array where the digest will be stored.
    //!
    void digest( uint8_t * d )
    {
      if ( m_finalized )
      {
        // If already finalized, just copy the state
        memcpy( d, m_state, static_cast<size_t>( digest_size() ) );
        return;
      }

      // Apply SHA3 padding
      m_apply_sha3_padding();

      // Squeeze
      memcpy( d, m_state, static_cast<size_t>( digest_size() ) );
      m_finalized = true;
    }

    //!
    //! \brief Returns a representation of the digest as a hexadecimal string.
    //!
    //! This method finalizes the hashing process and returns the resulting
    //! hash digest as a hexadecimal string. The caller takes ownership
    //! of the returned string.
    //!
    //! \return The hexadecimal string representation of the hash digest.
    //!
    std::string digest_in_hex()
    {
      // max digest_size = 100
      uint8_t bytes[100];

      // CHANGE: Uses bitshifting instead of multiplication
      digest( bytes );

      std::string hex;
      for ( int byte = 0; byte < digest_size(); byte++ )
      {
        // CHANGE: Uses bitshifting instead of multiplication
        hex += hexLookup[bytes[byte] >> 4];
        hex += hexLookup[bytes[byte] & 15];
      }
      return hex;
    }

    //!
    //! \brief Returns the size of the hash digest.
    //!
    //! This method returns the size of the hash output in bytes.
    //!
    //! \return The size of the digest in bytes.
    //!
    int digest_size() const { return m_digest_size; }

    //!
    //! \brief Reset the hash state for new computation.
    //!
    void reset() { m_reset(); }
  };
}  // namespace Utils

#endif
