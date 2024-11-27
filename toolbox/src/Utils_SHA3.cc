#include "Utils_SHA3.hh"
#include <algorithm>
#include <cstring>

using std::string;

// Circular rotate left
#define ROT_L( X, Y ) (( X << Y ) | ( X >> (64 - Y) ))
#define ROUNDS 24

/// For converting binary output to hexidecimal for printing
static char const * hexLookup = "0123456789abcdef";

static uint64_t const roundConstants[] = {
  0x0000000000000001,
  0x0000000000008082,
  0x800000000000808A,
  0x8000000080008000,
  0x000000000000808B,
  0x0000000080000001,
  0x8000000080008081,
  0x8000000000008009,
  0x000000000000008A,
  0x0000000000000088,
  0x0000000080008009,
  0x000000008000000A,
  0x000000008000808B,
  0x800000000000008B,
  0x8000000000008089,
  0x8000000000008003,
  0x8000000000008002,
  0x8000000000000080,
  0x000000000000800A,
  0x800000008000000A,
  0x8000000080008081,
  0x8000000000008080,
  0x0000000080000001,
  0x8000000080008008
};

namespace Utils {

  SHA3::SHA3( int digest_size ) : m_digest_size( digest_size ) {
    // zero the state
    // CHANGE: Now uses bit shifting instead of multiplication
    m_sponge_capacity = m_digest_size << 4; // x 16
    m_sponge_rate     = 1600 - m_sponge_capacity;
    m_reset();
  }

  ////////// Ingesting Data //////////

  void
  SHA3::hash( int const b ) {
    m_buffer_location[0] = static_cast<uint8_t>(b);
    m_buffer_location++;
    if ( m_buffer_location == reinterpret_cast<uint8_t*>(m_message_buffer_64) + (m_sponge_rate>>3) ) {
      m_buffer_location = reinterpret_cast<uint8_t*>(m_message_buffer_64);
      m_absorb_buffer();
    }
  }

  void
  SHA3::hash_string( const char * str ){
    int byte = 0;
    while( str[byte] != '\0' ){
      hash( static_cast<int>( static_cast<uint8_t>(str[byte]) ) );
      ++byte;
    }
  }

  void
  SHA3::hash_hex_string( const char *str ){
    int byte = 0;
    while( str[byte] != '\0' ){
      int f = str[byte];
      int s = str[byte+1];
      if      ( f >= 97 ) f -= 87; // lowercase
      else if ( f >= 65 ) f -= 55; // uppercase
      else                f -= 48; // numeric

      if      ( s >= 97 ) s -= 87; // lowercase
      else if ( s >= 65 ) s -= 55; // uppercase
      else                s -= 48; // numeric

      hash( (f << 4) | s );
      byte += 2;
    }
  }

  ////////// Expelling Data //////////

  void
  SHA3::digest( uint8_t * d ){
    // Pad with 10*1 padding
    m_buffer_location[0] = 1;
    m_buffer_location++;
    memset(
      m_buffer_location,
      0,
      size_t(reinterpret_cast<uint8_t*>(m_message_buffer_64) + (m_sponge_rate>>3) - m_buffer_location)
    );
    reinterpret_cast<uint8_t*>(m_message_buffer_64)[(m_sponge_rate >> 3) - 1] |= 0x80;
    m_absorb_buffer();

    // Squeeze
    memcpy( d, m_state, size_t(digest_size()) );
    m_reset(); // Ready the function to hash another message
  }

  string
  SHA3::digest_in_hex(){
    // max digest_size = 100
    // uint8_t *bytes = new uint8_t[ digest_size() ];
    // char *hex = new char[ (digest_size() << 1) + 1 ];

    uint8_t bytes[100];

    // CHANGE: Uses bitshifting instead of multiplication
    digest( bytes );

    string hex;
    for( int byte = 0; byte < digest_size(); byte++ ){
      // CHANGE: Uses bitshifting instead of multiplication
      hex += hexLookup[bytes[byte] >> 4];
      hex += hexLookup[bytes[byte] & 15];
    }
    return hex;
  }

  ////////// Internals //////////

  void
  SHA3::m_reset() {
    uint64_t * s = static_cast<uint64_t*>(m_state[0]);
    std::fill( s, s+25, 0 );
    m_buffer_location = reinterpret_cast<uint8_t*>(m_message_buffer_64);
  }

  void
  SHA3::m_absorb_buffer(){
    uint64_t *x = m_message_buffer_64;
    for( int i = 0; i*64 < m_sponge_rate; ++i ) m_state[i/5][i%5] ^= x[i]; // TODO: unroll
    m_perform_rounds( ROUNDS );
  }

  // CHANGE: Function changed to inline
  void
  SHA3::m_perform_rounds( int rounds ) {
    uint64_t b[5][5];
    uint64_t c[5];
    uint64_t d[5];

    for( int i = 0; i < rounds; i++ ) {

      //CHANGE: For loops change to pre-determined steps, reduces branching

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
      b[0][0] = m_state[0][0]; // rotate left by 0 bits
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
      m_state[0][0] = b[0][0] ^ ((~b[1][0]) & b[2][0]);
      m_state[1][0] = b[0][1] ^ ((~b[1][1]) & b[2][1]);
      m_state[2][0] = b[0][2] ^ ((~b[1][2]) & b[2][2]);
      m_state[3][0] = b[0][3] ^ ((~b[1][3]) & b[2][3]);
      m_state[4][0] = b[0][4] ^ ((~b[1][4]) & b[2][4]);

      m_state[0][1] = b[1][0] ^ ((~b[2][0]) & b[3][0]);
      m_state[1][1] = b[1][1] ^ ((~b[2][1]) & b[3][1]);
      m_state[2][1] = b[1][2] ^ ((~b[2][2]) & b[3][2]);
      m_state[3][1] = b[1][3] ^ ((~b[2][3]) & b[3][3]);
      m_state[4][1] = b[1][4] ^ ((~b[2][4]) & b[3][4]);

      m_state[0][2] = b[2][0] ^ ((~b[3][0]) & b[4][0]);
      m_state[1][2] = b[2][1] ^ ((~b[3][1]) & b[4][1]);
      m_state[2][2] = b[2][2] ^ ((~b[3][2]) & b[4][2]);
      m_state[3][2] = b[2][3] ^ ((~b[3][3]) & b[4][3]);
      m_state[4][2] = b[2][4] ^ ((~b[3][4]) & b[4][4]);

      m_state[0][3] = b[3][0] ^ ((~b[4][0]) & b[0][0]);
      m_state[1][3] = b[3][1] ^ ((~b[4][1]) & b[0][1]);
      m_state[2][3] = b[3][2] ^ ((~b[4][2]) & b[0][2]);
      m_state[3][3] = b[3][3] ^ ((~b[4][3]) & b[0][3]);
      m_state[4][3] = b[3][4] ^ ((~b[4][4]) & b[0][4]);

      m_state[0][4] = b[4][0] ^ ((~b[0][0]) & b[1][0]);
      m_state[1][4] = b[4][1] ^ ((~b[0][1]) & b[1][1]);
      m_state[2][4] = b[4][2] ^ ((~b[0][2]) & b[1][2]);
      m_state[3][4] = b[4][3] ^ ((~b[0][3]) & b[1][3]);
      m_state[4][4] = b[4][4] ^ ((~b[0][4]) & b[1][4]);

      // Iota step
      m_state[0][0] ^= roundConstants[i];
    }
  }

  ////////// Debugging Functions //////////

  void
  SHA3::m_print_message_buffer( ostream_type & stream ) const {
    uint8_t const * m_messageBuffer = reinterpret_cast<uint8_t const*>(m_message_buffer_64);
    stream << "mb = [ ";
    for( int i = 0; i < m_sponge_rate/8; ++i )
      stream << int(m_messageBuffer[i]) << ' ';
    stream << "]\n";
  }

  void
  SHA3::m_print_sponge( ostream_type & stream ) const {
    stream << "s = [ " << std::hex;
    for( int i = 0; i < 5; ++i )
      for( int j = 0; j < 5; ++j )
        stream << m_state[i][j] << ' ';
    stream << std::dec << "]\n";
  }
}
