//
//  main.cpp
//  genconjson
//
//  Created by Nicola Dal Bianco on 08/11/17.
//  Copyright © 2017 Nicola Dal Bianco. All rights reserved.
//

#ifndef GC_JSON_HANDLER_HH
#define GC_JSON_HANDLER_HH

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#include "GenericContainerJson.hh"

#ifdef USE_SYSTEM_JSON
  #include <rapidjson/reader.h>
#else
  #include "rapidjson/reader.h"
#endif

#endif

namespace GC_namespace {

  class GenericContainerJsonHandler : public rapidjson::BaseReaderHandler< rapidjson::UTF8<>, GenericContainerJsonHandler > {

  public:

    GenericContainerJsonHandler() = delete;

    GenericContainerJsonHandler ( GenericContainer & gc_output, GenericContainer const & gc_options );

    // BaseReaderHandler implementation
    bool Null();
    bool Bool ( bool b );
    bool Int ( int i );
    bool Uint ( unsigned int u );
    bool Int64 ( int64_t i );
    bool Uint64 ( uint64_t u );
    bool Double ( double d );
    bool String ( const char * str, rapidjson::SizeType length, bool copy );
    bool StartObject();
    bool Key ( const char * str, rapidjson::SizeType length, bool copy );
    bool EndObject ( rapidjson::SizeType member_count );
    bool StartArray();
    bool EndArray ( rapidjson::SizeType member_count );

  protected:

    typedef std::pair<GenericContainer *, bool> stack_entry;

    std::vector<stack_entry> _gc_stack = {};

    GenericContainer * getCurrentGCPointer() const;

    GCJsonMatrixOrder _mat_order;

    std::string _im_unit;

    void setCurrentGCPointerArrayType( bool is_array );
    bool isCurrentGCPointerArrayType() const;
    void parseStringObjectToCurrentGc( std::string const & str );
    void advanceCurrentGCPointer();
    void finalizeArrayProcess();

  };

}

#endif
