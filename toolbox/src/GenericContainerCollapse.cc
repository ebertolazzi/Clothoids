/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: GenericContainerCollapse.cc
//

#include "GenericContainer/GenericContainer.hh"
#include <iomanip>
#include <cmath>
#include <ctgmath>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

namespace GC_namespace {

  /*
  //   _                                 _
  //  | |_ ___     _   _  __ _ _ __ ___ | |
  //  | __/ _ \   | | | |/ _` | '_ ` _ \| |
  //  | || (_) |  | |_| | (_| | | | | | | |
  //   \__\___/____\__, |\__,_|_| |_| |_|_|
  //         |_____|___/
  */

  void
  GenericContainer::collapse() {
    switch (m_data_type) {
      case GC_type::MAP: {
        map_type & M{*m_data.m};
        for ( auto & m : M ) m.second.collapse();
        break;
      }
      case GC_type::VECTOR: {
        // posso collassare solo i vettori
        vector_type & v{*m_data.v};
        TypeAllowed max_sub_type{GC_type::NOTYPE};
        bool can_collapse{true};
        int_type last_nelem{-1};
        for ( auto & vi : v ) {
          vi.collapse();
          // get type
          TypeAllowed ti{vi.get_type()};
          if ( max_sub_type < ti ) max_sub_type = ti;
          int_type ne = int_type(vi.get_num_elements());
          if ( last_nelem >= 0 && can_collapse ) {
            can_collapse = last_nelem == ne;
          }
          last_nelem = ne;
        }
        if ( !can_collapse ) break;
        switch ( max_sub_type ) {
          case GC_type::BOOL: {
            vec_bool_type tmp; tmp.reserve(v.size());
            for ( auto & vi : v ) tmp.push_back( vi.get_bool() );
            *this = tmp; // sovrascrive
            break;
          }
          case GC_type::INTEGER: { vec_int_type tmp;     this->copyto_vec_int(tmp);     *this = tmp; break; }
          case GC_type::LONG:    { vec_long_type tmp;    this->copyto_vec_long(tmp);    *this = tmp; break; }
          case GC_type::REAL:    { vec_real_type tmp;    this->copyto_vec_real(tmp);    *this = tmp; break; }
          case GC_type::COMPLEX: { vec_complex_type tmp; this->copyto_vec_complex(tmp); *this = tmp; break; }
          case GC_type::VEC_INTEGER: {
            unsigned NR{unsigned(last_nelem)};
            unsigned NC{unsigned(v.size())};
            mat_int_type M(NR,NC);
            for ( unsigned j{0}; j < NC; ++j ) {
              vec_int_type C; v[j].copyto_vec_int(C);
              for ( unsigned i{0}; i < NR; ++i ) M(i,j) = C[i];
            }
            *this = M;
            break;
          }
          case GC_type::VEC_LONG: {
            unsigned NR{unsigned(last_nelem)};
            unsigned NC{unsigned(v.size())};
            mat_long_type M(NR,NC);
            for ( unsigned j{0}; j < NC; ++j ) {
              vec_long_type C; v[j].copyto_vec_long(C);
              for ( unsigned i{0}; i < NR; ++i ) M(i,j) = C[i];
            }
            *this = M;
            break;
          }
          case GC_type::VEC_REAL: {
            unsigned NR{unsigned(last_nelem)};
            unsigned NC{unsigned(v.size())};
            mat_real_type M(NR,NC);
            for ( unsigned j{0}; j < NC; ++j ) {
              vec_real_type C; v[j].copyto_vec_real(C);
              for ( unsigned i{0}; i < NR; ++i ) M(i,j) = C[i];
            }
            *this = M;
            break;
          }
          case GC_type::VEC_COMPLEX: {
            unsigned NR{unsigned(last_nelem)};
            unsigned NC{unsigned(v.size())};
            mat_complex_type M(NR,NC);
            for ( unsigned j{0}; j < NC; ++j ) {
              vec_complex_type C; v[j].copyto_vec_complex(C);
              for ( unsigned i{0}; i < NR; ++i ) M(i,j) = C[i];
            }
            *this = M;
            break;
          }
          default:
          break;
        }
        break;
      }
      default:
        break;
    }
  }
}

//
// eof: GenericContainerCollapse.cc
//
