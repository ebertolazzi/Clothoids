#pragma once

#ifndef UTILS_SHA3_HH
#define UTILS_SHA3_HH

#include "Utils.hh"

#include <string>
#include <iostream>
#include <cstdint>

namespace Utils {

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
  class SHA3 {

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

    //!
    //! \brief Constructor for the SHA3 class
    //!
    //! Initializes the SHA3 hashing algorithm with the specified digest size.
    //!
    //! \param digest_size The size of the desired hash output in bytes (e.g., 224, 256, 384, or 512).
    //!
    SHA3( int digest_size );

    //! \brief Destructor for the SHA3 class.
    ~SHA3() { }

    //!
    //! \brief Adds an entire string to the message.
    //!
    //! This method appends the given string to the internal message buffer
    //! for hashing. It processes the string as a sequence of bytes.
    //!
    //! \param str The null-terminated string of bytes to add to the hash.
    //!
    void hash_string( char const * str );

    //!
    //! \brief Adds an entire hexadecimal string to the message.
    //!
    //! This method appends the given hexadecimal string to the internal
    //! message buffer. Each pair of hex digits is converted to a byte.
    //!
    //! \param str The null-terminated hexadecimal string of bytes to add to the hash.
    //!
    void hash_hex_string( char const * str );

    //!
    //! \brief Returns a representation of the digest as a hexadecimal string.
    //!
    //! This method finalizes the hashing process and returns the resulting
    //! hash digest as a hexadecimal string. The caller takes ownership
    //! of the returned string.
    //!
    //! \return The hexadecimal string representation of the hash digest.
    //!
    string digest_in_hex();


    //!
    //! \brief Returns the size of the hash digest.
    //!
    //! This method returns the size of the hash output in bytes.
    //!
    //! \return The size of the digest in bytes.
    //!
    int digest_size() const { return m_digest_size; }

    //!
    //! \brief Hash a specified number of bytes.
    //!
    //! This method processes the provided byte input to update the hash state.
    //!
    //! \param b The number of bytes to hash.
    //!
    void hash( const int b );

    //!
    //! \brief Retrieve the digest and store it in the provided byte array.
    //!
    //! This method finalizes the hashing process and stores the resulting
    //! hash digest in the provided byte array.
    //!
    //! \param d Pointer to the byte array where the digest will be stored.
    //!
    void digest( uint8_t * d );
  };
}

#endif
