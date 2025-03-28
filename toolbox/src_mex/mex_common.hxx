// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
BaseCurve *
get_base_pointer( string const & kind, mxArray const * arg ) {
  // mexPrintf( "get_base_pointer -> for: %s\n", kind.c_str() );
  if ( kind == "linesegment" )
    return Utils::mex_convert_mx_to_ptr<LineSegment>(arg);
  else if ( kind == "circlearc" )
    return Utils::mex_convert_mx_to_ptr<CircleArc>(arg);
  else if ( kind == "biarc" )
    return Utils::mex_convert_mx_to_ptr<Biarc>(arg);
  else if ( kind == "clothoidcurve" )
    return Utils::mex_convert_mx_to_ptr<ClothoidCurve>(arg);
  else if ( kind == "polyline" )
    return Utils::mex_convert_mx_to_ptr<PolyLine>(arg);
  else if ( kind == "biarclist"  )
    return Utils::mex_convert_mx_to_ptr<BiarcList>(arg);
  else if ( kind == "clothoidlist"  )
    return Utils::mex_convert_mx_to_ptr<ClothoidList>(arg);
  else if ( kind == "dubins"  )
    return Utils::mex_convert_mx_to_ptr<Dubins>(arg);

  UTILS_MEX_ASSERT( false, "in get_base_pointer type '{}' unknown\n", kind );
  return nullptr;
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
bool
do_is_ISO( mxArray const * plhs, char const msg[] ) {
  UTILS_MEX_ASSERT0( mxIsChar(plhs), msg );
  string cmd = mxArrayToString( plhs );
  return cmd == "ISO";
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_length(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('length',OBJ[,offs,'ISO'/'SAE']): "
  UTILS_MEX_ASSERT(
    nrhs == 2 || nrhs == 3 || nrhs == 4,
    CMD "expected 2, 3 or 4 inputs, nrhs = {}\n", nrhs
  );
  UTILS_MEX_ASSERT(
    nlhs == 1,
    CMD "expected 1 output, nlhs = {}\n", nlhs
  );

  if ( nrhs == 2 ) {
    Utils::mex_set_scalar_value( arg_out_0, ptr->length() );
  } else {
    real_type offs;
    offs = Utils::mex_get_scalar_value(
      arg_in_2,
      CMD "`offs` expected to be a real scalar"
    );
    bool ISO = true;
    if ( nrhs == 4 ) ISO = do_is_ISO( arg_in_3, CMD " last argument must be a string");
    real_type len = ISO ? ptr->length_ISO( offs ) : ptr->length_SAE( offs );
    Utils::mex_set_scalar_value( arg_out_0, len );
  }

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_make_a_copy(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  #define CMD CMD_BASE "('make_a_copy',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  G2LIB_CLASS const * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);
  arg_out_0 = Utils::mex_convert_ptr_to_mx<G2LIB_CLASS>(new G2LIB_CLASS( *ptr ));
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_build2(
  int nlhs, mxArray       *[],
  int nrhs, mxArray const *prhs[]
) {

  #define CMD CMD_BASE "('build2',OBJ,OBJ1,type1): "
  UTILS_MEX_ASSERT( nrhs == 4, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );

  UTILS_MEX_ASSERT0( mxIsChar(arg_in_3), CMD "'type' argument must be a string" );
  string kind = mxArrayToString( arg_in_3 );
  Utils::to_lower(kind);

  G2LIB_CLASS     * ptr0 = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);
  BaseCurve const * ptr1 = get_base_pointer( kind, arg_in_2 );

  if ( ptr1->type() == CurveType::LINE )
    ptr0->build( *static_cast<LineSegment const *>(ptr1) );
  else if ( ptr1->type() == CurveType::CIRCLE )
    ptr0->build( *static_cast<CircleArc const *>(ptr1) );
  else if ( ptr1->type() == CurveType::BIARC )
    ptr0->build( *static_cast<Biarc const *>(ptr1) );
  else if ( ptr1->type() == CurveType::CLOTHOID )
    ptr0->build( *static_cast<ClothoidCurve const *>(ptr1) );
  else if ( ptr1->type() == CurveType::POLYLINE )
    ptr0->build( *static_cast<PolyLine const *>(ptr1) );
  else if ( ptr1->type() == CurveType::BIARC_LIST )
    ptr0->build( *static_cast<BiarcList const *>(ptr1) );
  else if ( ptr1->type() == CurveType::CLOTHOID_LIST )
    ptr0->build( *static_cast<ClothoidList const *>(ptr1) );
  else if ( ptr1->type() == CurveType::DUBINS )
    ptr0->build( *static_cast<Dubins const *>(ptr1) );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_delete(
  int nlhs, mxArray       *[],
  int nrhs, mxArray const *prhs[]
) {
  #define CMD CMD_BASE "('delete',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );
  // Destroy the C++ object
  Utils::mex_destroy_object<G2LIB_CLASS>( arg_in_1 );
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_bbox(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  #define CMD CMD_BASE "('bbox',OBJ[,offs,'ISO'/'SAE']): "
  UTILS_MEX_ASSERT(
    nrhs >= 2 && nrhs <= 4,
    CMD "expected 2, 3 or 4 inputs, nrhs = {}\n", nrhs
  );
  UTILS_MEX_ASSERT(
    nlhs == 4,
    CMD "expected 4 output, nlhs = {}\n", nlhs
  );

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  real_type xmin, ymin, xmax, ymax;
  if ( nrhs >= 3 ) {
    real_type offs = Utils::mex_get_scalar_value(
      arg_in_2, CMD "`offs` expected to be a real scalar"
    );
    bool ISO = true;
    if ( nrhs == 4 ) ISO = do_is_ISO( arg_in_3, CMD " last argument must be a string");
    if ( ISO ) ptr->bbox_ISO( offs, xmin, ymin, xmax, ymax );
    else       ptr->bbox_SAE( offs, xmin, ymin, xmax, ymax );
  } else {
    ptr->bbox( xmin, ymin, xmax, ymax );
  }

  Utils::mex_set_scalar_value( arg_out_0, xmin );
  Utils::mex_set_scalar_value( arg_out_1, ymin );
  Utils::mex_set_scalar_value( arg_out_2, xmax );
  Utils::mex_set_scalar_value( arg_out_3, ymax );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_bbTriangles( int nlhs, mxArray       * plhs[],
                int nrhs, mxArray const * prhs[] ) {

  #define CMD CMD_BASE "('bbTriangles',OBJ[,max_angle,max_size,offs,'ISO'/'SAE']): "

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  UTILS_MEX_ASSERT(
    nrhs >= 2 && nrhs <= 6,
    CMD "expected 2 up to 6 inputs, nrhs = {}\n", nrhs
  );
  UTILS_MEX_ASSERT(
    nlhs == 3,
    CMD "expected 3 output, nlhs = {}\n", nlhs
  );

  real_type max_angle = Utils::m_pi/18;
  real_type max_size  = 1e100;
  real_type offs      = 0;
  if ( nrhs >= 3 )
    max_angle = Utils::mex_get_scalar_value(
      arg_in_2, CMD "`max_angle` expected to be a real scalar"
    );
  if ( nrhs >= 4 )
    max_size = Utils::mex_get_scalar_value(
      arg_in_3, CMD "`max_size` expected to be a real scalar"
    );
  if ( nrhs >= 5 )
    offs = Utils::mex_get_scalar_value(
      arg_in_4, CMD "`offs` expected to be a real scalar"
    );

  bool ISO = true;
  if ( nrhs == 6 ) ISO = do_is_ISO( arg_in_5, CMD " last argument must be a string");

  std::vector<Triangle2D> tvec;
  if ( nrhs >= 5 ) {
    if ( ISO ) {
      ptr->bb_triangles_ISO( offs, tvec, max_angle, max_size, 0 );
    } else {
      ptr->bb_triangles_SAE( offs, tvec, max_angle, max_size, 0 );
    }
  } else {
    ptr->bb_triangles( tvec, max_angle, max_size, 0 );
  }

  mwSize nt = tvec.size();

  double * p0 = Utils::mex_create_matrix_value( arg_out_0, 2, nt );
  double * p1 = Utils::mex_create_matrix_value( arg_out_1, 2, nt );
  double * p2 = Utils::mex_create_matrix_value( arg_out_2, 2, nt );

  for ( mwSize i = 0; i < nt; ++i ) {
    Triangle2D const & t = tvec[i];
    *p0++ = t.x1();
    *p0++ = t.y1();
    *p1++ = t.x2();
    *p1++ = t.y2();
    *p2++ = t.x3();
    *p2++ = t.y3();
  }

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_change_origin(
  int nlhs, mxArray       *[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('change_origin',OBJ,x0,y0): "
  UTILS_MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );

  real_type new_x0, new_y0;
  new_x0 = Utils::mex_get_scalar_value( arg_in_2, CMD "`x0` expected to be a real scalar" );
  new_y0 = Utils::mex_get_scalar_value( arg_in_3, CMD "`y0` expected to be a real scalar" );

  ptr->change_origin( new_x0, new_y0 );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_translate(
  int nlhs, mxArray       *[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('translate',OBJ,t0,t0): "
  UTILS_MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );

  real_type tx, ty;
  tx = Utils::mex_get_scalar_value( arg_in_2, CMD "`tx` expected to be a real scalar" );
  ty = Utils::mex_get_scalar_value( arg_in_3, CMD "`ty` expected to be a real scalar" );

  ptr->translate( tx, ty );
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_rotate(
  int nlhs, mxArray       *[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('rotate',OBJ,angle,cx,cy): "
  UTILS_MEX_ASSERT( nrhs == 5, CMD "expected 5 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );

  real_type angle, cx, cy;
  angle = Utils::mex_get_scalar_value( arg_in_2, CMD "`angle` expected to be a real scalar" );
  cx    = Utils::mex_get_scalar_value( arg_in_3, CMD "`cx` expected to be a real scalar" );
  cy    = Utils::mex_get_scalar_value( arg_in_4, CMD "`cy` expected to be a real scalar" );

  ptr->rotate( angle, cx, cy );
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_scale(
  int nlhs, mxArray       *[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('scale',OBJ,sc): "
  UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );
  real_type sc = Utils::mex_get_scalar_value( arg_in_2, CMD "`sc` expected to be a real scalar" );
  ptr->scale(sc);
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_reverse(
  int nlhs, mxArray       *[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('reverse',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );

  ptr->reverse();
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_trim(
  int nlhs, mxArray       *[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('trim',OBJ,s_begin,s_end): "
  UTILS_MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );

  real_type s_begin, s_end;
  s_begin = Utils::mex_get_scalar_value( arg_in_2, CMD "`s_begin` expected to be a real scalar" );
  s_end   = Utils::mex_get_scalar_value( arg_in_3, CMD "`s_end` expected to be a real scalar" );

  ptr->trim( s_begin, s_end );
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_closest_point(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('closest_point',OBJ,qx,qy[,offs,'ISO'/'SAE']): "

  UTILS_MEX_ASSERT(
    nrhs >= 4 || nrhs <= 6,
    CMD "expected 4, 5 or 6 input, nrhs = {}\n", nrhs
  );
  UTILS_MEX_ASSERT(
    nlhs == 1 || nlhs == 6,
    CMD "expected 1 pr 6 outputs, nlhs = {}\n", nlhs
  );

  mwSize nrx, ncx, nry, ncy;
  real_type const * qx;
  real_type const * qy;
  qx = Utils::mex_matrix_pointer(
    arg_in_2, nrx, ncx,
    CMD "`qx` expected to be a real vector/matrix"
  );
  qy = Utils::mex_matrix_pointer(
    arg_in_3, nry, ncy,
    CMD "`qy` expected to be a real vector/matrix"
  );

  UTILS_MEX_ASSERT(
    nrx == nry && ncx == ncy,
    CMD "`qx` and `qy` expected to be of the same size, found\n"
    "size(qx) = {} x {} size(qy) = {} x {}\n",
    nrx, ncx, nry, ncy
  );

  real_type * x     = nullptr;
  real_type * y     = nullptr;
  real_type * s     = nullptr;
  real_type * t     = nullptr;
  int32_t   * iflag = nullptr;
  real_type * dst   = nullptr;

  if ( nlhs == 1 ) {

    mxArray * mx_x;
    mxArray * mx_y;
    mxArray * mx_s;
    mxArray * mx_t;
    mxArray * mx_iflag;
    mxArray * mx_dst;

    x      = Utils::mex_create_matrix_value( mx_x,     nrx, ncx );
    y      = Utils::mex_create_matrix_value( mx_y,     nrx, ncx );
    s      = Utils::mex_create_matrix_value( mx_s,     nrx, ncx );
    t      = Utils::mex_create_matrix_value( mx_t,     nrx, ncx );
    iflag  = Utils::mex_create_matrix_int32( mx_iflag, nrx, ncx );
    dst    = Utils::mex_create_matrix_value( mx_dst,   nrx, ncx );

    static char const * fieldnames[] = {
      "x", "y", "s", "t", "iflag", "dst"
    };

    arg_out_0 = mxCreateStructMatrix(1,1,6,fieldnames);
    mxSetFieldByNumber( arg_out_0, 0, 0, mx_x     );
    mxSetFieldByNumber( arg_out_0, 0, 1, mx_y     );
    mxSetFieldByNumber( arg_out_0, 0, 2, mx_s     );
    mxSetFieldByNumber( arg_out_0, 0, 3, mx_t     );
    mxSetFieldByNumber( arg_out_0, 0, 4, mx_iflag );
    mxSetFieldByNumber( arg_out_0, 0, 5, mx_dst   );

  } else {
    x     = Utils::mex_create_matrix_value( arg_out_0, nrx, ncx );
    y     = Utils::mex_create_matrix_value( arg_out_1, nrx, ncx );
    s     = Utils::mex_create_matrix_value( arg_out_2, nrx, ncx );
    t     = Utils::mex_create_matrix_value( arg_out_3, nrx, ncx );
    iflag = Utils::mex_create_matrix_int32( arg_out_4, nrx, ncx );
    dst   = Utils::mex_create_matrix_value( arg_out_5, nrx, ncx );
  }

  mwSize size = nrx * ncx;
  if ( nrhs >= 5 ) {
    real_type offs = Utils::mex_get_scalar_value(
      arg_in_4, CMD "`offs` expected to be a real scalar"
    );
    bool ISO = true;
    if ( nrhs == 6 ) ISO = do_is_ISO( arg_in_5, CMD " last argument must be a string");
    if ( ISO ) {
      for ( mwSize i = 0; i < size; ++i )
        *iflag++ = ptr->closest_point_ISO(
          *qx++,  *qy++, offs, *x++, *y++, *s++, *t++, *dst++
        );
    } else {
      for ( mwSize i = 0; i < size; ++i )
        *iflag++ = ptr->closest_point_SAE(
          *qx++,  *qy++, offs, *x++, *y++, *s++, *t++, *dst++
        );
    }
  } else {
    for ( mwSize i = 0; i < size; ++i )
      *iflag++ = ptr->closest_point_ISO(
        *qx++, *qy++, *x++, *y++, *s++, *t++, *dst++
      );
  }

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_distance(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('distance',OBJ,qx,qy[,offs,'ISO'/'SAE']): "

  UTILS_MEX_ASSERT(
    nrhs >= 4 || nrhs <= 6,
    CMD "expected 4, 5 or 6 input, nrhs = {}\n", nrhs
  );
  UTILS_MEX_ASSERT(
    nlhs == 1,
    CMD "expected 1 output, nlhs = {}\n", nlhs
  );

  mwSize nrx, ncx, nry, ncy;
  real_type const * qx;
  real_type const * qy;
  qx = Utils::mex_matrix_pointer(
    arg_in_2, nrx, ncx,
    CMD "`qx` expected to be a real vector/matrix"
  );
  qy = Utils::mex_matrix_pointer(
    arg_in_3, nry, ncy,
    CMD "`qy` expected to be a real vector/matrix"
  );

  UTILS_MEX_ASSERT(
    nrx == nry && ncx == ncy,
    CMD "`qx` and `qy` expected to be of the same size, found\n"
    "size(qx) = {} x {} size(qy) = {} x {}\n",
    nrx, ncx, nry, ncy
  );

  real_type * dst = Utils::mex_create_matrix_value( arg_out_0, nrx, ncx );

  mwSize size = nrx * ncx;
  if ( nrhs >= 5 ) {
    real_type offs = Utils::mex_get_scalar_value(
      arg_in_4, CMD "`offs` expected to be a real scalar"
    );
    bool ISO = true;
    if ( nrhs == 6 ) ISO = do_is_ISO( arg_in_5, CMD " last argument must be a string");
    if ( ISO ) {
      for ( mwSize i = 0; i < size; ++i )
        *dst++ = ptr->distance_ISO( *qx++, *qy++, offs );
    } else {
      for ( mwSize i = 0; i < size; ++i )
        *dst++ = ptr->distance_SAE( *qx++, *qy++, offs );
    }
  } else {
    for ( mwSize i = 0; i < size; ++i )
      *dst++ = ptr->distance( *qx++, *qy++ );
  }

  #undef CMD
}

static
void
do_collision(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  #define CMD CMD_BASE "('collision',OBJ,OBJ1,type[,offs,offs1,'ISO'/'SAE']): "
  UTILS_MEX_ASSERT( nrhs >= 4 || nrhs <= 7, CMD "expected 4 to 7 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  UTILS_MEX_ASSERT0( mxIsChar(arg_in_3), CMD "'type' argument must be a string" );
  string kind = mxArrayToString( arg_in_3 );
  Utils::to_lower(kind);

  G2LIB_CLASS * ptr0 = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);
  BaseCurve   * ptr1 = get_base_pointer( kind, arg_in_2 );

  G2LIB_DEBUG_MESSAGE(
    "do_collision ADDRS: {}, {}\n", fmt::ptr(ptr0), fmt::ptr(ptr1)
  );

  if ( nrhs == 4 ) {
    Utils::mex_set_scalar_bool( arg_out_0, G2lib::collision( ptr0, ptr1 ) );
  } else {
    real_type offs     = Utils::mex_get_scalar_value( arg_in_4, CMD "`offs` expected to be a real scalar" );
    real_type offs_obj = Utils::mex_get_scalar_value( arg_in_5, CMD "`offs_obj` expected to be a real scalar" );
    bool ISO = true;
    if ( nrhs == 7 ) ISO = do_is_ISO( arg_in_6, CMD " last argument must be a string");
    bool ok;
    if ( ISO ) ok = G2lib::collision_ISO( ptr0, offs, ptr1, offs_obj );
    else       ok = G2lib::collision_SAE( ptr0, offs, ptr1, offs_obj );
    Utils::mex_set_scalar_bool( arg_out_0, ok );
  }

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_intersect(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  #define CMD CMD_BASE "( 'intersect', OBJ, OBJ1, type, [,offs, offs1, 'ISO'/'SAE'] ): "
  UTILS_MEX_ASSERT( nrhs >= 4 || nrhs <= 7, CMD "expected 4 to 7 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = {}\n", nlhs );

  UTILS_MEX_ASSERT0( mxIsChar(arg_in_3), CMD "'type' argument must be a string" );
  string kind = mxArrayToString( arg_in_3 );
  Utils::to_lower(kind);

  G2LIB_CLASS * ptr0 = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);
  BaseCurve   * ptr1 = get_base_pointer( kind, arg_in_2 );

  G2LIB_DEBUG_MESSAGE(
    "do_intersect ADDRS: {}, {}\n", fmt::ptr(ptr0), fmt::ptr(ptr1)
  );

  IntersectList ilist;
  if ( nrhs == 4 ) {
    intersect( ptr0, ptr1, ilist );
  } else {
    real_type offs, offs_obj;
    offs     = Utils::mex_get_scalar_value( arg_in_4, CMD "`offs` expected to be a real scalar" );
    offs_obj = Utils::mex_get_scalar_value( arg_in_5, CMD "`offs_obj` expected to be a real scalar" );
    bool ISO = true;
    if ( nrhs == 7 ) ISO = do_is_ISO( arg_in_6, CMD " last argument must be a string" );
    if ( ISO ) intersect_ISO( ptr0, offs, ptr1, offs_obj, ilist );
    else       intersect_SAE( ptr0, offs, ptr1, offs_obj, ilist );
  }

  double * pS1 = Utils::mex_create_matrix_value( arg_out_0, ilist.size(), 1 );
  double * pS2 = Utils::mex_create_matrix_value( arg_out_1, ilist.size(), 1 );
  for ( mwSize i=0; i < ilist.size(); ++i ) {
    Ipair const & ip = ilist[i];
    *pS1++ = ip.first;
    *pS2++ = ip.second;
  }

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_findST(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('findST',OBJ,x,y[,'ISO'/'SAE']): "
  UTILS_MEX_ASSERT( nrhs == 4 || nrhs == 5, CMD "expected 4 or 5 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = {}\n", nlhs );
  mwSize nrx, ncx, nry, ncy;
  real_type const * x;
  real_type const * y;
  x = Utils::mex_matrix_pointer(
    arg_in_2, nrx, ncx,
    CMD "`x` expected to be a real vector/matrix"
  );
  y = Utils::mex_matrix_pointer(
    arg_in_3, nry, ncy,
    CMD "`y` expected to be a real vector/matrix"
  );
  UTILS_MEX_ASSERT(
    nrx == nry && ncx == ncy,
    CMD "`x` and `y` expected to be of the same size, found\n"
    "size(x) = {} x {} size(y) = {} x {}\n",
    nrx, ncx, nry, ncy
  );

  real_type * s = Utils::mex_create_matrix_value( arg_out_0, nrx, ncx );
  real_type * t = Utils::mex_create_matrix_value( arg_out_1, nrx, ncx );

  bool ISO = true;
  if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string");

  mwSize size = nrx*ncx;
  if ( ISO ) {
    for ( mwSize i = 0; i < size; ++i )
      ptr->findST_ISO( *x++, *y++, *s++, *t++ );
  } else {
    for ( mwSize i = 0; i < size; ++i )
      ptr->findST_SAE( *x++, *y++, *s++, *t++ );
  }
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_info(
  int nlhs, mxArray       *[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('info',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO outputs, nlhs = {}\n", nlhs );
  ptr->info( std::cout );
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_eval(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  if ( nrhs >= 4 ) {

    #define CMD CMD_BASE "('eval',OBJ,s,[offs,'ISO'/'SAE'): "

    mwSize size, sizet;
    real_type const * s;
    real_type const * t;
    s = Utils::mex_vector_pointer( arg_in_2, size,  CMD "`s` expected to be a real vector" );
    t = Utils::mex_vector_pointer( arg_in_3, sizet, CMD "`t` expected to be a real vector" );

    bool ISO = true;
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string" );

    UTILS_MEX_ASSERT(
      size == sizet || size == 1 || sizet ==1,
      CMD " size(s) = {} must be equal to size(t) = {} or size(s|t) == 1\n",
      size, sizet
    );

    mwSize incs = size  == 1 ? 0 : 1;
    mwSize inct = sizet == 1 ? 0 : 1;
    mwSize npts = size > sizet ? size : sizet;

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = Utils::mex_create_matrix_value( arg_out_0, 2, size );
      if ( ISO ) {
        LOOPXY1 ptr->eval_ISO( *s, *t, pXY[0], pXY[1] );
      } else {
        LOOPXY1 ptr->eval_SAE( *s, *t, pXY[0], pXY[1] );
      }
    } else if ( nlhs == 2 ) {
      real_type *pX = Utils::mex_create_matrix_value( arg_out_0, 1, size );
      real_type *pY = Utils::mex_create_matrix_value( arg_out_1, 1, size );
      if ( ISO ) {
        LOOPXY2 ptr->eval_ISO( *s, *t, *pX, *pY );
      } else {
        LOOPXY2 ptr->eval_SAE( *s, *t, *pX, *pY );
      }
    } else {
      UTILS_MEX_ASSERT(
        nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = {}\n", nlhs
      );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD

  } else if ( nrhs == 3 ) {

    #define CMD CMD_BASE "('eval',OBJ,s): "

    mwSize npts;
    real_type const * s;
    s = Utils::mex_vector_pointer( arg_in_2, npts, CMD "`s` expected to be a real vector" );

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, ++s, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, ++s, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = Utils::mex_create_matrix_value( arg_out_0, 2, npts );
      LOOPXY1 ptr->eval( *s, pXY[0], pXY[1] );
    } else if ( nlhs == 2 ) {
      real_type *pX = Utils::mex_create_matrix_value( arg_out_0, 1, npts );
      real_type *pY = Utils::mex_create_matrix_value( arg_out_1, 1, npts );
      LOOPXY2 ptr->eval( *s, *pX, *pY );
    } else {
      UTILS_MEX_ASSERT(
        nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = {}\n", nlhs
      );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD
  }
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_eval_D(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  if ( nrhs == 4 || nrhs == 5 ) {

    #define CMD CMD_BASE "('eval_D',OBJ,s,t[,'ISO'/'SAE']): "

    mwSize size, sizet;
    real_type const * s;
    real_type const * t;
    s = Utils::mex_vector_pointer( arg_in_2, size, CMD "`s` expected to be a real vector" );
    t = Utils::mex_vector_pointer( arg_in_3, sizet, CMD "`t` expected to be a real vector" );

    bool ISO = true;
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string");

    UTILS_MEX_ASSERT(
      size == sizet || size == 1 || sizet ==1,
      CMD " size(s) = {} must be equal to size(t) = {} or size(s|t) == 1\n",
      size, sizet
    );

    mwSize incs = size  == 1 ? 0 : 1;
    mwSize inct = sizet == 1 ? 0 : 1;
    mwSize npts = size > sizet ? size : sizet;

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = Utils::mex_create_matrix_value( arg_out_0, 2,size );
      if ( ISO ) {
        LOOPXY1 ptr->eval_ISO_D( *s, *t, pXY[0], pXY[1] );
      } else {
        LOOPXY1 ptr->eval_SAE_D( *s, *t, pXY[0], pXY[1] );
      }
    } else if ( nlhs == 2 ) {
      real_type *pX = Utils::mex_create_matrix_value( arg_out_0, 1,size );
      real_type *pY = Utils::mex_create_matrix_value( arg_out_1, 1,size );
      if ( ISO ) {
        LOOPXY2 ptr->eval_ISO_D( *s, *t, *pX, *pY );
      } else {
        LOOPXY2 ptr->eval_SAE_D( *s, *t, *pX, *pY );
      }
    } else {
      UTILS_MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = {}\n", nlhs );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD

  } else if ( nrhs == 3 ) {

    #define CMD CMD_BASE "('eval_D',OBJ,s): "

    mwSize npts;
    real_type const * s;
    s = Utils::mex_vector_pointer(
      arg_in_2, npts, CMD "`s` expected to be a real vector"
    );

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, ++s, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, ++s, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = Utils::mex_create_matrix_value( arg_out_0, 2, npts );
      LOOPXY1 ptr->eval_D( *s, pXY[0], pXY[1] );
    } else if ( nlhs == 2 ) {
      real_type *pX = Utils::mex_create_matrix_value( arg_out_0, 1, npts );
      real_type *pY = Utils::mex_create_matrix_value( arg_out_1, 1, npts );
      LOOPXY2 ptr->eval_D( *s, *pX, *pY );
    } else {
      UTILS_MEX_ASSERT(
        nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = {}\n", nlhs
      );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD
  }
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_eval_DD(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  if ( nrhs == 4 || nrhs == 5 ) {

    #define CMD CMD_BASE "('eval_DD',OBJ,s,t[,'ISO'/'SAE']): "

    mwSize size, sizet;
    real_type const * s;
    real_type const * t;
    s = Utils::mex_vector_pointer( arg_in_2, size, CMD "`s` expected to be a real vector" );
    t = Utils::mex_vector_pointer( arg_in_3, sizet, CMD "`t` expected to be a real vector" );

    UTILS_MEX_ASSERT(
      size == sizet || size == 1 || sizet ==1,
      CMD " size(s) = {} must be equal to size(t) = {} or size(s|t) == 1\n",
      size, sizet
    );

    mwSize incs = size  == 1 ? 0 : 1;
    mwSize inct = sizet == 1 ? 0 : 1;
    mwSize npts = size > sizet ? size : sizet;

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, ++pX, ++pY )

    bool ISO = true;
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string");

    if ( nlhs == 1 ) {
      real_type *pXY = Utils::mex_create_matrix_value( arg_out_0, 2,size );
      if ( ISO ) {
        LOOPXY1 ptr->eval_ISO_DD( *s, *t, pXY[0], pXY[1] );
      } else {
        LOOPXY1 ptr->eval_SAE_DD( *s, *t, pXY[0], pXY[1] );
      }
    } else if ( nlhs == 2 ) {
      real_type *pX = Utils::mex_create_matrix_value( arg_out_0, 1,size );
      real_type *pY = Utils::mex_create_matrix_value( arg_out_1, 1,size );
      if ( ISO ) {
        LOOPXY2 ptr->eval_ISO_DD( *s, *t, *pX, *pY );
      } else {
        LOOPXY2 ptr->eval_SAE_DD( *s, *t, *pX, *pY );
      }
    } else {
      UTILS_MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = {}\n", nlhs );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD

  } else if ( nrhs == 3 ) {

    #define CMD CMD_BASE "('eval_DD',OBJ,s): "

    mwSize npts;
    real_type const * s;
    s = Utils::mex_vector_pointer(
      arg_in_2, npts, CMD "`s` expected to be a real vector"
    );

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, ++s, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, ++s, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = Utils::mex_create_matrix_value( arg_out_0, 2, npts );
      LOOPXY1 ptr->eval_DD( *s, pXY[0], pXY[1] );
    } else if ( nlhs == 2 ) {
      real_type *pX = Utils::mex_create_matrix_value( arg_out_0, 1, npts );
      real_type *pY = Utils::mex_create_matrix_value( arg_out_1, 1, npts );
      LOOPXY2 ptr->eval_DD( *s, *pX, *pY );
    } else {
      UTILS_MEX_ASSERT(
        nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = {}\n", nlhs
      );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD
  }
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_eval_DDD(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  if ( nrhs == 4 || nrhs == 5 ) {

    #define CMD CMD_BASE "('eval_DDD',OBJ,s,t[,'ISO'/'SAE']): "

    mwSize size, sizet;
    real_type const * s;
    real_type const * t;
    s = Utils::mex_vector_pointer( arg_in_2, size, CMD "`s` expected to be a real vector" );
    t = Utils::mex_vector_pointer( arg_in_3, sizet, CMD "`t` expected to be a real vector" );

    UTILS_MEX_ASSERT(
      size == sizet || size == 1 || sizet ==1,
      CMD " size(s) = {} must be equal to size(t) = {} or size(s|t) == 1\n",
      size, sizet
    );

    mwSize incs = size  == 1 ? 0 : 1;
    mwSize inct = sizet == 1 ? 0 : 1;
    mwSize npts = size > sizet ? size : sizet;

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, ++pX, ++pY )

    bool ISO = true;
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string");

    if ( nlhs == 1 ) {
      real_type *pXY = Utils::mex_create_matrix_value( arg_out_0, 2,size );
      if ( ISO ) {
        LOOPXY1 ptr->eval_ISO_DDD( *s, *t, pXY[0], pXY[1] );
      } else {
        LOOPXY1 ptr->eval_SAE_DDD( *s, *t, pXY[0], pXY[1] );
      }
    } else if ( nlhs == 2 ) {
      real_type *pX = Utils::mex_create_matrix_value( arg_out_0, 1,size );
      real_type *pY = Utils::mex_create_matrix_value( arg_out_1, 1,size );
      if ( ISO ) {
        LOOPXY2 ptr->eval_ISO_DDD( *s, *t, *pX, *pY );
      } else {
        LOOPXY2 ptr->eval_SAE_DDD( *s, *t, *pX, *pY );
      }
    } else {
      UTILS_MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = {}\n", nlhs );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD

  } else if ( nrhs == 3 ) {

    #define CMD CMD_BASE "('eval_DDD',OBJ,s): "

    mwSize npts;
    real_type const * s;
    s = Utils::mex_vector_pointer( arg_in_2, npts, CMD "`s` expected to be a real vector" );

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, ++s, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, ++s, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = Utils::mex_create_matrix_value( arg_out_0, 2, npts );
      LOOPXY1 ptr->eval_DDD( *s, pXY[0], pXY[1] );
    } else if ( nlhs == 2 ) {
      real_type *pX = Utils::mex_create_matrix_value( arg_out_0, 1, npts );
      real_type *pY = Utils::mex_create_matrix_value( arg_out_1, 1, npts );
      LOOPXY2 ptr->eval_DDD( *s, *pX, *pY );
    } else {
      UTILS_MEX_ASSERT(
        nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = {}\n", nlhs
      );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD
  }

}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_evaluate(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  if ( nrhs >= 4 && nrhs <= 5) {

    #define CMD CMD_BASE "('evaluate',OBJ,s,[offs,'ISO'/'SAE']): "

    mwSize size, sizet;
    real_type const * s;
    real_type const * t;
    s = Utils::mex_vector_pointer( arg_in_2, size, CMD "`s` expected to be a real vector" );
    t = Utils::mex_vector_pointer( arg_in_3, sizet, CMD "`t` expected to be a real vector" );
    bool ISO = true;
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string");

    UTILS_MEX_ASSERT(
      size == sizet || size == 1 || sizet ==1,
      CMD " size(s) = {} must be equal to size(t) = {} or size(s|t) == 1\n",
      size, sizet
    );

    mwSize incs = size  == 1 ? 0 : 1;
    mwSize inct = sizet == 1 ? 0 : 1;
    mwSize npts = size > sizet ? size : sizet;

    #define LOOPXY1 \
      for ( mwSize i = 0; i < npts;  \
            ++i, s += incs, t += inct, pXY += 4 )
    #define LOOPXY2 \
      for ( mwSize i = 0; i < npts; \
            ++i, s += incs, t += inct, ++pTH, ++pK, ++pX, ++pY  )

    if ( nlhs == 1 ) {
      real_type *pXY = Utils::mex_create_matrix_value( arg_out_0, 4, size );
      if ( ISO ) {
        LOOPXY1 ptr->evaluate_ISO( *s, *t, pXY[0], pXY[1], pXY[2], pXY[3] );
      } else {
        LOOPXY1 ptr->evaluate_SAE( *s, *t, pXY[0], pXY[1], pXY[2], pXY[3] );
      }
    } else if ( nlhs == 4 ) {
      real_type *pX  = Utils::mex_create_matrix_value( arg_out_0, 1, size );
      real_type *pY  = Utils::mex_create_matrix_value( arg_out_1, 1, size );
      real_type *pTH = Utils::mex_create_matrix_value( arg_out_2, 1, size );
      real_type *pK  = Utils::mex_create_matrix_value( arg_out_3, 1, size );
      if ( ISO ) {
        LOOPXY2 ptr->evaluate_ISO( *s, *t, *pTH, *pK, *pX, *pY );
      } else {
        LOOPXY2 ptr->evaluate_SAE( *s, *t, *pTH, *pK, *pX, *pY );
      }
    } else {
      UTILS_MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = {}\n", nlhs );
    }
    #undef LOOPXY1
    #undef LOOPXY2

  } else if ( nrhs == 3 ) {

    mwSize npts;
    real_type const * s;
    s = Utils::mex_vector_pointer( arg_in_2, npts, CMD "`s` expected to be a real vector" );

    #define LOOPXY1 \
      for ( mwSize i = 0; i < npts;  \
            ++i, ++s, pXY += 4 )
    #define LOOPXY2 \
      for ( mwSize i = 0; i < npts; \
            ++i, ++s, ++pTH, ++pK, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = Utils::mex_create_matrix_value( arg_out_0, 4, npts );
      LOOPXY1 ptr->evaluate( *s, pXY[0], pXY[1], pXY[2], pXY[3] );
    } else if ( nlhs == 4 ) {
      real_type *pX  = Utils::mex_create_matrix_value( arg_out_0, 1, npts );
      real_type *pY  = Utils::mex_create_matrix_value( arg_out_1, 1, npts );
      real_type *pTH = Utils::mex_create_matrix_value( arg_out_2, 1, npts );
      real_type *pK  = Utils::mex_create_matrix_value( arg_out_3, 1, npts );
      LOOPXY2 ptr->evaluate( *s, *pTH, *pK, *pX, *pY );
    } else {
      UTILS_MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = {}\n", nlhs );
    }

    #undef LOOPXY1
    #undef LOOPXY2
  }
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_theta(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('theta',OBJ,s): "
  UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  mwSize size;
  real_type const * s;
  s = Utils::mex_vector_pointer( arg_in_2, size, CMD "`s` expected to be a real vector" );

  real_type *theta = Utils::mex_create_matrix_value( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i )  *theta++ = ptr->theta( *s++ );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_theta_D(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('theta_D',OBJ,s): "
  UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  mwSize size;
  real_type const * s;
  s = Utils::mex_vector_pointer( arg_in_2, size, CMD "`s` expected to be a real vector" );

  real_type *theta = Utils::mex_create_matrix_value( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i ) *theta++ = ptr->theta_D( *s++ );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_theta_DD(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('theta_DD',OBJ,s): "
  UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  mwSize size;
  real_type const * s;
  s = Utils::mex_vector_pointer( arg_in_2, size, CMD "`s` expected to be a real vector" );

  real_type *theta = Utils::mex_create_matrix_value( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i )
    *theta++ = ptr->theta_DD( *s++ );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_theta_DDD(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('theta_DDD',OBJ,s): "
  UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  mwSize size;
  real_type const * s;
  s = Utils::mex_vector_pointer( arg_in_2, size, CMD "`s` expected to be a real vector" );

  real_type *theta = Utils::mex_create_matrix_value( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i )  *theta++ = ptr->theta_DDD( *s++ );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_kappa(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('kappa',OBJ,s): "
  UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  mwSize size;
  real_type const * s;
  s = Utils::mex_vector_pointer( arg_in_2, size, CMD "`s` expected to be a real vector" );

  real_type *kappa = Utils::mex_create_matrix_value( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i ) *kappa++ = ptr->kappa( *s++ );
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_kappa_D(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('kappa_D',OBJ,s): "
  UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  mwSize size;
  real_type const * s;
  s = Utils::mex_vector_pointer(
    arg_in_2, size, CMD "`s` expected to be a real vector"
  );

  real_type *kappa_D = Utils::mex_create_matrix_value( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i ) *kappa_D++ = ptr->kappa_D( *s++ );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_kappa_DD(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('kappa_DD',OBJ,s): "
  UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  mwSize size;
  real_type const * s;
  s = Utils::mex_vector_pointer(
    arg_in_2, size, CMD "`s` expected to be a real vector"
  );

  real_type *kappa_DD = Utils::mex_create_matrix_value( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i )
    *kappa_DD++ = ptr->kappa_DD( *s++ );

  #undef CMD

}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_xy_begin(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('xy_begin',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = {}\n", nlhs );

  Utils::mex_set_scalar_value( arg_out_0, ptr->x_begin() );
  Utils::mex_set_scalar_value( arg_out_1, ptr->y_begin() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_x_begin(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('x_begin',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  Utils::mex_set_scalar_value( arg_out_0, ptr->x_begin() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_y_begin(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('y_begin',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  Utils::mex_set_scalar_value( arg_out_0, ptr->y_begin() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_theta_begin(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('theta_begin',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  Utils::mex_set_scalar_value( arg_out_0, ptr->theta_begin() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_kappa_begin(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('kappa_begin',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  Utils::mex_set_scalar_value( arg_out_0, ptr->kappa_begin() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_xy_end(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('xy_end',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = {}\n", nlhs );

  Utils::mex_set_scalar_value( arg_out_0, ptr->x_end() );
  Utils::mex_set_scalar_value( arg_out_1, ptr->y_end() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_x_end(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('x_end',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  Utils::mex_set_scalar_value( arg_out_0, ptr->x_end() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_y_end(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('y_end',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  Utils::mex_set_scalar_value( arg_out_0, ptr->y_end() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_theta_end(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('theta_end',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  Utils::mex_set_scalar_value( arg_out_0, ptr->theta_end() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_kappa_end(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = Utils::mex_convert_mx_to_ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('kappa_end',OBJ): "
  UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

  Utils::mex_set_scalar_value( arg_out_0, ptr->kappa_end() );

  #undef CMD
}

static
void
do_yesAABBtree(
  int nlhs, mxArray       *[],
  int nrhs, mxArray const *[]
) {

  #define CMD CMD_BASE "('yesAABBtree'): "
  UTILS_MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );

  G2lib::yesAABBtree();

  #undef CMD
}

static
void
do_noAABBtree(
  int nlhs, mxArray       *[],
  int nrhs, mxArray const *[]
) {

  #define CMD CMD_BASE "('noAABBtree'): "
  UTILS_MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = {}\n", nrhs );
  UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );

  G2lib::noAABBtree();

  #undef CMD
}

#define CMD_MAP_FUN                 \
{"length",do_length},               \
{"delete",do_delete},               \
{"make_a_copy",do_make_a_copy},     \
{"build2",do_build2},               \
{"bbox",do_bbox},                   \
{"bbTriangles",do_bbTriangles},     \
{"change_origin",do_change_origin}, \
{"translate",do_translate},         \
{"rotate",do_rotate},               \
{"scale",do_scale},                 \
{"reverse",do_reverse},             \
{"trim",do_trim},                   \
{"distance",do_distance},           \
{"closest_point",do_closest_point}, \
{"collision",do_collision},         \
{"intersect",do_intersect},         \
{"findST",do_findST},               \
{"info",do_info},                   \
{"evaluate",do_evaluate},           \
{"eval",do_eval},                   \
{"eval_D",do_eval_D},               \
{"eval_DD",do_eval_DD},             \
{"eval_DDD",do_eval_DDD},           \
{"theta",do_theta},                 \
{"theta_D",do_theta_D},             \
{"theta_DD",do_theta_DD},           \
{"theta_DDD",do_theta_DDD},         \
{"kappa",do_kappa},                 \
{"kappa_D",do_kappa_D},             \
{"kappa_DD",do_kappa_DD},           \
{"xy_begin",do_xy_begin},           \
{"x_begin",do_x_begin},             \
{"y_begin",do_y_begin},             \
{"theta_begin",do_theta_begin},     \
{"kappa_begin",do_kappa_begin},     \
{"xy_end",do_xy_end},               \
{"x_end",do_x_end},                 \
{"y_end",do_y_end},                 \
{"theta_end",do_theta_end},         \
{"kappa_end",do_kappa_end},         \
{"yesAABBtree",do_yesAABBtree},     \
{"noAABBtree",do_noAABBtree}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
