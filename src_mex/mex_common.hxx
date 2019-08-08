// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
bool
do_is_ISO( mxArray const * plhs, char const msg[] ) {
  MEX_ASSERT( mxIsChar(plhs), msg );
  string cmd = mxArrayToString( plhs );
  return cmd == "ISO";
}

static
void
do_length(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('length',OBJ[,offs,'ISO'/'SAE']): "
  MEX_ASSERT(
    nrhs == 2 || nrhs == 3 || nrhs == 4,
    CMD "expected 2, 3 or 4 inputs, nrhs = " << nrhs
  );
  MEX_ASSERT(
    nlhs == 1,
    CMD "expected 1 output, nlhs = " << nlhs
  );

  if ( nrhs == 2 ) {
    setScalarValue( arg_out_0, ptr->length() );
  } else {
    real_type offs;
    offs = getScalarValue(
      arg_in_2,
      CMD "`offs` expected to be a real scalar"
    );
    bool ISO = true;
    if ( nrhs == 4 ) ISO = do_is_ISO( arg_in_3, CMD " last argument must be a string");
    real_type len = ISO ? ptr->length_ISO( offs ) : ptr->length_SAE( offs );
    setScalarValue( arg_out_0, len );
  }

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_copy(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('copy',OBJ,OBJ1): "
  MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

  G2LIB_CLASS const * CC = convertMat2Ptr<G2LIB_CLASS>(arg_in_2);
  ptr->copy(*CC);

  #undef CMD
}

static
void
do_delete(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define CMD CMD_BASE "('delete',OBJ): "
  MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
  MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
  // Destroy the C++ object
  destroyObject<G2LIB_CLASS>( arg_in_1 );
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
  MEX_ASSERT(
    nrhs >= 2 && nrhs <= 4,
    CMD "expected 2, 3 or 4 inputs, nrhs = " << nrhs
  );
  MEX_ASSERT(
    nlhs == 4,
    CMD "expected 4 output, nlhs = " << nlhs
  );

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  real_type xmin, ymin, xmax, ymax;
  if ( nrhs >= 3 ) {
    real_type offs = getScalarValue(
      arg_in_2, CMD "`offs` expected to be a real scalar"
    );
    bool ISO = true;
    if ( nrhs == 4 ) ISO = do_is_ISO( arg_in_3, CMD " last argument must be a string");
    if ( ISO ) ptr->bbox_ISO( offs, xmin, ymin, xmax, ymax );
    else       ptr->bbox_SAE( offs, xmin, ymin, xmax, ymax );
  } else {
    ptr->bbox( xmin, ymin, xmax, ymax );
  }

  setScalarValue( arg_out_0, xmin );
  setScalarValue( arg_out_1, ymin );
  setScalarValue( arg_out_2, xmax );
  setScalarValue( arg_out_3, ymax );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_change_origin(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('changeOrigin',OBJ,x0,y0): "
  MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

  real_type new_x0, new_y0;
  new_x0 = getScalarValue(
    arg_in_2, CMD "`x0` expected to be a real scalar"
  );
  new_y0 = getScalarValue(
    arg_in_3, CMD "`y0` expected to be a real scalar"
  );

  ptr->changeOrigin( new_x0, new_y0 );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_translate(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('translate',OBJ,t0,t0): "
  MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

  real_type tx, ty;
  tx = getScalarValue(
    arg_in_2, CMD "`tx` expected to be a real scalar"
  );
  ty = getScalarValue(
    arg_in_3, CMD "`ty` expected to be a real scalar"
  );

  ptr->translate( tx, ty );
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_rotate(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('rotate',OBJ,angle,cx,cy): "
  MEX_ASSERT( nrhs == 5, CMD "expected 5 inputs, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

  real_type angle, cx, cy;
  angle = getScalarValue( arg_in_2, CMD "`angle` expected to be a real scalar" );
  cx    = getScalarValue( arg_in_3, CMD "`cx` expected to be a real scalar" );
  cy    = getScalarValue( arg_in_4, CMD "`cy` expected to be a real scalar" );

  ptr->rotate( angle, cx, cy );
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_scale(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('scale',OBJ,sc): "
  MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
  real_type sc = getScalarValue( arg_in_2, CMD "`sc` expected to be a real scalar" );
  ptr->scale(sc);
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_reverse(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('reverse',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

  ptr->reverse();
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_trim(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('trim',OBJ,s_begin,s_end): "
  MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

  real_type s_begin, s_end;
  s_begin = getScalarValue(
    arg_in_2, CMD "`s_begin` expected to be a real scalar"
  );
  s_end = getScalarValue(
    arg_in_3, CMD "`s_end` expected to be a real scalar"
  );

  ptr->trim( s_begin, s_end );
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_closestPoint(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('closestPoint',OBJ,qx,qy[,offs,'ISO'/'SAE']): "

  MEX_ASSERT(
    nrhs >= 4 || nrhs <= 6,
    CMD "expected 4, 5 or 6 input, nrhs = " << nrhs
  );
  MEX_ASSERT(
    nlhs == 6,
    CMD "expected 6 output, nlhs = " << nlhs
  );

  mwSize nrx, ncx, nry, ncy;
  real_type const * qx;
  real_type const * qy;
  qx = getMatrixPointer(
    arg_in_2, nrx, ncx,
    CMD "`qx` expected to be a real vector/matrix"
  );
  qy = getMatrixPointer(
    arg_in_3, nry, ncy,
    CMD "`qy` expected to be a real vector/matrix"
  );

  MEX_ASSERT(
    nrx == nry && ncx == ncy,
    CMD "`qx` and `qy` expected to be of the same size, found size(qx) = " <<
    nrx << " x " << nry << " size(qy) = " << nry << " x " << ncy
  );

  real_type * x     = createMatrixValue( arg_out_0, nrx, ncx );
  real_type * y     = createMatrixValue( arg_out_1, nrx, ncx );
  real_type * s     = createMatrixValue( arg_out_2, nrx, ncx );
  real_type * t     = createMatrixValue( arg_out_3, nrx, ncx );
  int32_t   * iflag = createMatrixInt32( arg_out_4, nrx, ncx );
  real_type * dst   = createMatrixValue( arg_out_5, nrx, ncx );

  mwSize size = nrx * ncx;
  if ( nrhs >= 5 ) {
    real_type offs = getScalarValue(
      arg_in_4, CMD "`offs` expected to be a real scalar"
    );
    bool ISO = true;
    if ( nrhs == 6 ) ISO = do_is_ISO( arg_in_5, CMD " last argument must be a string");
    if ( ISO ) {
      for ( mwSize i = 0; i < size; ++i )
        *iflag++ = ptr->closestPoint_ISO(
          *qx++,  *qy++, offs, *x++, *y++, *s++, *t++, *dst++
        );
    } else {
      for ( mwSize i = 0; i < size; ++i )
        *iflag++ = ptr->closestPoint_SAE(
          *qx++,  *qy++, offs, *x++, *y++, *s++, *t++, *dst++
        );
    }
  } else {
    for ( mwSize i = 0; i < size; ++i )
      *iflag++ = ptr->closestPoint_ISO(
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

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('distance',OBJ,qx,qy[,offs,'ISO'/'SAE']): "

  MEX_ASSERT(
    nrhs >= 4 || nrhs <= 6,
    CMD "expected 4, 5 or 6 input, nrhs = " << nrhs
  );
  MEX_ASSERT(
    nlhs == 1,
    CMD "expected 1 output, nlhs = " << nlhs
  );

  mwSize nrx, ncx, nry, ncy;
  real_type const * qx;
  real_type const * qy;
  qx = getMatrixPointer(
    arg_in_2, nrx, ncx,
    CMD "`qx` expected to be a real vector/matrix"
  );
  qy = getMatrixPointer(
    arg_in_3, nry, ncy,
    CMD "`qy` expected to be a real vector/matrix"
  );

  MEX_ASSERT(
    nrx == nry && ncx == ncy,
    CMD "`qx` and `qy` expected to be of the same size, found size(qx) = " <<
    nrx << " x " << nry << " size(qy) = " << nry << " x " << ncy
  );

  real_type * dst = createMatrixValue( arg_out_0, nrx, ncx );

  mwSize size = nrx * ncx;
  if ( nrhs >= 5 ) {
    real_type offs = getScalarValue(
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

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_collision(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  #define CMD CMD_BASE "('collision',OBJ,OBJ1,type[,offs,offs1,'ISO'/'SAE']): "
  MEX_ASSERT(
    nrhs >= 4 || nrhs <= 7,
    CMD "expected 4 to 7 inputs, nrhs = " << nrhs
  );
  MEX_ASSERT(
    nlhs == 1,
    CMD "expected 1 output, nlhs = " << nlhs
  );

  MEX_ASSERT( mxIsChar(arg_in_3), CMD "'type' argument must be a string" );
  string kind = mxArrayToString( arg_in_3 );

  G2LIB_CLASS * ptr  = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);
  BaseCurve   * ptr1 = nullptr;

  if      ( kind == "LineSegment"   ) ptr1 = convertMat2Ptr<LineSegment>(arg_in_2);
  else if ( kind == "CircleArc"     ) ptr1 = convertMat2Ptr<CircleArc>(arg_in_2);
  else if ( kind == "BiArc"         ) ptr1 = convertMat2Ptr<Biarc>(arg_in_2);
  else if ( kind == "ClothoidCurve" ) ptr1 = convertMat2Ptr<ClothoidCurve>(arg_in_2);
  else if ( kind == "ClothoidList"  ) ptr1 = convertMat2Ptr<ClothoidList>(arg_in_2);
  else if ( kind == "PolyLine"      ) ptr1 = convertMat2Ptr<PolyLine>(arg_in_2);
  else {
    MEX_ASSERT( false, CMD "'type '" << "' unsupported" );
  }

  if ( nrhs == 4 ) {
    setScalarBool( arg_out_0, collision( *ptr, *ptr1 ) );
  } else {
    real_type offs, offs_obj;
    offs = getScalarValue(
      arg_in_4, CMD "`offs` expected to be a real scalar"
    );
    offs_obj = getScalarValue(
      arg_in_5, CMD "`offs_obj` expected to be a real scalar"
    );
    bool ISO = true;
    if ( nrhs == 7 ) ISO = do_is_ISO( arg_in_6, CMD " last argument must be a string");
    bool ok;
    if ( ISO ) ok = collision_ISO( *ptr, offs, *ptr1, offs_obj );
    else       ok = collision_SAE( *ptr, offs, *ptr1, offs_obj );
    setScalarBool( arg_out_0, ok );
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

  #define CMD CMD_BASE "('intersect',OBJ,OBJ1,type,[,offs,offs1,'ISO'/'SAE']): "
  MEX_ASSERT(
    nrhs >= 4 || nrhs <= 7,
    CMD "expected 4 to 7 inputs, nrhs = " << nrhs
  );
  MEX_ASSERT(
    nlhs == 2,
    CMD "expected 2 output, nlhs = " << nlhs
  );

  MEX_ASSERT( mxIsChar(arg_in_3), CMD "'type' argument must be a string" );
  string kind = mxArrayToString( arg_in_3 );

  G2LIB_CLASS * ptr  = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);
  BaseCurve   * ptr1 = nullptr;

  if      ( kind == "LineSegment"   ) ptr1 = convertMat2Ptr<LineSegment>(arg_in_2);
  else if ( kind == "CircleArc"     ) ptr1 = convertMat2Ptr<CircleArc>(arg_in_2);
  else if ( kind == "BiArc"         ) ptr1 = convertMat2Ptr<Biarc>(arg_in_2);
  else if ( kind == "ClothoidCurve" ) ptr1 = convertMat2Ptr<ClothoidCurve>(arg_in_2);
  else if ( kind == "ClothoidList"  ) ptr1 = convertMat2Ptr<ClothoidList>(arg_in_2);
  else if ( kind == "PolyLine"      ) ptr1 = convertMat2Ptr<PolyLine>(arg_in_2);
  else {
    MEX_ASSERT( false, CMD "'type '" << "' unsupported" );
  }

  IntersectList ilist;
  if ( nrhs == 4 ) {
    intersect( *ptr, *ptr1, ilist, false );
  } else {
    real_type offs, offs_obj;
    offs = getScalarValue(
      arg_in_4, CMD "`offs` expected to be a real scalar"
    );
    offs_obj = getScalarValue(
      arg_in_5, CMD "`offs_obj` expected to be a real scalar"
    );
    bool ISO = true;
    if ( nrhs == 7 ) ISO = do_is_ISO( arg_in_6, CMD " last argument must be a string");
    if ( ISO ) intersect_ISO( *ptr, offs, *ptr1, offs_obj, ilist, false );
    else       intersect_SAE( *ptr, offs, *ptr1, offs_obj, ilist, false );
  }

  double * pS1 = createMatrixValue( arg_out_0, ilist.size(), 1 );
  double * pS2 = createMatrixValue( arg_out_1, ilist.size(), 1 );
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

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('findST',OBJ,x,y[,'ISO'/'SAE']): "
  MEX_ASSERT( nrhs == 4 || nrhs == 5, CMD "expected 4 or 5 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = " << nlhs );
  mwSize nrx, ncx, nry, ncy;
  real_type const * x;
  real_type const * y;
  x = getMatrixPointer(
    arg_in_2, nrx, ncx,
    CMD "`x` expected to be a real vector/matrix"
  );
  y = getMatrixPointer(
    arg_in_3, nry, ncy,
    CMD "`y` expected to be a real vector/matrix"
  );
  MEX_ASSERT( 
    nrx == nry && ncx == ncy,
    CMD "`x` and `y` expected to be of the same size, found size(x) = " <<
    nrx << " x " << nry << " size(y) = " << nry << " x " << ncy
  );

  real_type * s = createMatrixValue( arg_out_0, nrx, ncx );
  real_type * t = createMatrixValue( arg_out_1, nrx, ncx );

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
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('info',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 0, CMD "expected NO outputs, nlhs = " << nlhs );
  ptr->info(cout);
  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_eval(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  if ( nrhs >= 4 ) {

    #define CMD CMD_BASE "('eval',OBJ,s,[offs,'ISO'/'SAE'): "

    mwSize size, sizet;
    real_type const * s;
    real_type const * t;
    s = getVectorPointer(
      arg_in_2, size, CMD "`s` expected to be a real vector"
    );
    t = getVectorPointer(
      arg_in_3, sizet, CMD "`t` expected to be a real vector"
    );

    bool ISO = true;
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string");

    MEX_ASSERT(
      size == sizet || size == 1 || sizet ==1,
      CMD " size(s) = " << size << " must be equal to size(t) = " << 
      sizet << " or size(s|t) == 1"
    );

    mwSize incs = size  == 1 ? 0 : 1;
    mwSize inct = sizet == 1 ? 0 : 1;
    mwSize npts = max(size,sizet);

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = createMatrixValue( arg_out_0, 2,size );
      if ( ISO ) {
        LOOPXY1 ptr->eval_ISO( *s, *t, pXY[0], pXY[1] );
      } else {
        LOOPXY1 ptr->eval_SAE( *s, *t, pXY[0], pXY[1] );
      }
    } else if ( nlhs == 2 ) {
      real_type *pX = createMatrixValue( arg_out_0, 1,size );
      real_type *pY = createMatrixValue( arg_out_1, 1,size );
      if ( ISO ) {
        LOOPXY2 ptr->eval_ISO( *s, *t, *pX, *pY );
      } else {
        LOOPXY2 ptr->eval_SAE( *s, *t, *pX, *pY );
      }
    } else {
      MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD

  } else if ( nrhs == 3 ) {

    #define CMD CMD_BASE "('eval',OBJ,s): "

    mwSize npts;
    real_type const * s;
    s = getVectorPointer(
      arg_in_2, npts,
      CMD "`s` expected to be a real vector"
    );

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, ++s, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, ++s, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = createMatrixValue( arg_out_0, 2, npts );
      LOOPXY1 ptr->eval( *s, pXY[0], pXY[1] );
    } else if ( nlhs == 2 ) {
      real_type *pX = createMatrixValue( arg_out_0, 1, npts );
      real_type *pY = createMatrixValue( arg_out_1, 1, npts );
      LOOPXY2 ptr->eval( *s, *pX, *pY );
    } else {
      MEX_ASSERT(
        nlhs == 0,
        CMD "expected 1 or 2 outputs, nlhs = " << nlhs
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

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  if ( nrhs == 4 || nrhs == 5 ) {

    #define CMD CMD_BASE "('eval_D',OBJ,s,t[,'ISO'/'SAE']): "

    mwSize size, sizet;
    real_type const * s;
    real_type const * t;
    s = getVectorPointer(
      arg_in_2, size, CMD "`s` expected to be a real vector"
    );
    t = getVectorPointer(
      arg_in_3, sizet, CMD "`t` expected to be a real vector"
    );

    bool ISO = true;
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string");

    MEX_ASSERT(
      size == sizet || size == 1 || sizet ==1,
      CMD " size(s) = " << size <<
      " must be equal to size(t) = " << sizet <<
      " or size(s|t) == 1"
    );

    mwSize incs = size  == 1 ? 0 : 1;
    mwSize inct = sizet == 1 ? 0 : 1;
    mwSize npts = max(size,sizet);

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = createMatrixValue( arg_out_0, 2,size );
      if ( ISO ) {
        LOOPXY1 ptr->eval_ISO_D( *s, *t, pXY[0], pXY[1] );
      } else {
        LOOPXY1 ptr->eval_SAE_D( *s, *t, pXY[0], pXY[1] );
      }
    } else if ( nlhs == 2 ) {
      real_type *pX = createMatrixValue( arg_out_0, 1,size );
      real_type *pY = createMatrixValue( arg_out_1, 1,size );
      if ( ISO ) {
        LOOPXY2 ptr->eval_ISO_D( *s, *t, *pX, *pY );
      } else {
        LOOPXY2 ptr->eval_SAE_D( *s, *t, *pX, *pY );
      }
    } else {
      MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD

  } else if ( nrhs == 3 ) {

    #define CMD CMD_BASE "('eval_D',OBJ,s): "

    mwSize npts;
    real_type const * s;
    s = getVectorPointer(
      arg_in_2, npts, CMD "`s` expected to be a real vector"
    );

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, ++s, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, ++s, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = createMatrixValue( arg_out_0, 2, npts );
      LOOPXY1 ptr->eval_D( *s, pXY[0], pXY[1] );
    } else if ( nlhs == 2 ) {
      real_type *pX = createMatrixValue( arg_out_0, 1, npts );
      real_type *pY = createMatrixValue( arg_out_1, 1, npts );
      LOOPXY2 ptr->eval_D( *s, *pX, *pY );
    } else {
      MEX_ASSERT(
        nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs
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

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  if ( nrhs == 4 || nrhs == 5 ) {

    #define CMD CMD_BASE "('eval_DD',OBJ,s,t[,'ISO'/'SAE']): "

    mwSize size, sizet;
    real_type const * s;
    real_type const * t;
    s = getVectorPointer(
      arg_in_2, size, CMD "`s` expected to be a real vector"
    );
    t = getVectorPointer(
      arg_in_3, sizet, CMD "`t` expected to be a real vector" 
    );

    MEX_ASSERT(
      size == sizet || size == 1 || sizet ==1,
      CMD " size(s) = " << size <<
      " must be equal to size(t) = " << sizet <<
      " or size(s|t) == 1"
    );

    mwSize incs = size  == 1 ? 0 : 1;
    mwSize inct = sizet == 1 ? 0 : 1;
    mwSize npts = max(size,sizet);

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, ++pX, ++pY )

    bool ISO = true;
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string");

    if ( nlhs == 1 ) {
      real_type *pXY = createMatrixValue( arg_out_0, 2,size );
      if ( ISO ) {
        LOOPXY1 ptr->eval_ISO_DD( *s, *t, pXY[0], pXY[1] );
      } else {
        LOOPXY1 ptr->eval_SAE_DD( *s, *t, pXY[0], pXY[1] );
      }
    } else if ( nlhs == 2 ) {
      real_type *pX = createMatrixValue( arg_out_0, 1,size );
      real_type *pY = createMatrixValue( arg_out_1, 1,size );
      if ( ISO ) {
        LOOPXY2 ptr->eval_ISO_DD( *s, *t, *pX, *pY );
      } else {
        LOOPXY2 ptr->eval_SAE_DD( *s, *t, *pX, *pY );
      }
    } else {
      MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD

  } else if ( nrhs == 3 ) {

    #define CMD CMD_BASE "('eval_DD',OBJ,s): "

    mwSize npts;
    real_type const * s;
    s = getVectorPointer(
      arg_in_2, npts, CMD "`s` expected to be a real vector"
    );

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, ++s, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, ++s, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = createMatrixValue( arg_out_0, 2, npts );
      LOOPXY1 ptr->eval_DD( *s, pXY[0], pXY[1] );
    } else if ( nlhs == 2 ) {
      real_type *pX = createMatrixValue( arg_out_0, 1, npts );
      real_type *pY = createMatrixValue( arg_out_1, 1, npts );
      LOOPXY2 ptr->eval_DD( *s, *pX, *pY );
    } else {
      MEX_ASSERT(
        nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs
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

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  if ( nrhs == 4 || nrhs == 5 ) {

    #define CMD CMD_BASE "('eval_DDD',OBJ,s,t[,'ISO'/'SAE']): "

    mwSize size, sizet;
    real_type const * s;
    real_type const * t;
    s = getVectorPointer(
      arg_in_2, size, CMD "`s` expected to be a real vector"
    );
    t = getVectorPointer(
      arg_in_3, sizet, CMD "`t` expected to be a real vector"
    );

    MEX_ASSERT(
      size == sizet || size == 1 || sizet ==1,
      CMD " size(s) = " << size <<
      " must be equal to size(t) = " << sizet <<
      " or size(s|t) == 1"
    );

    mwSize incs = size  == 1 ? 0 : 1;
    mwSize inct = sizet == 1 ? 0 : 1;
    mwSize npts = max(size,sizet);

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, ++pX, ++pY )

    bool ISO = true;
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string");

    if ( nlhs == 1 ) {
      real_type *pXY = createMatrixValue( arg_out_0, 2,size );
      if ( ISO ) {
        LOOPXY1 ptr->eval_ISO_DDD( *s, *t, pXY[0], pXY[1] );
      } else {
        LOOPXY1 ptr->eval_SAE_DDD( *s, *t, pXY[0], pXY[1] );
      }
    } else if ( nlhs == 2 ) {
      real_type *pX = createMatrixValue( arg_out_0, 1,size );
      real_type *pY = createMatrixValue( arg_out_1, 1,size );
      if ( ISO ) {
        LOOPXY2 ptr->eval_ISO_DDD( *s, *t, *pX, *pY );
      } else {
        LOOPXY2 ptr->eval_SAE_DDD( *s, *t, *pX, *pY );
      }
    } else {
      MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD

  } else if ( nrhs == 3 ) {

    #define CMD CMD_BASE "('eval_DDD',OBJ,s): "

    mwSize npts;
    real_type const * s;
    s = getVectorPointer(
      arg_in_2, npts, CMD "`s` expected to be a real vector"
    );

    #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, ++s, pXY += 2 )
    #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, ++s, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = createMatrixValue( arg_out_0, 2, npts );
      LOOPXY1 ptr->eval_DDD( *s, pXY[0], pXY[1] );
    } else if ( nlhs == 2 ) {
      real_type *pX = createMatrixValue( arg_out_0, 1, npts );
      real_type *pY = createMatrixValue( arg_out_1, 1, npts );
      LOOPXY2 ptr->eval_DDD( *s, *pX, *pY );
    } else {
      MEX_ASSERT(
        nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs
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

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  if ( nrhs >= 4 && nrhs <= 5) {

    #define CMD CMD_BASE "('evaluate',OBJ,s,offs[,'ISO'/'SAE']): "

    mwSize size, sizet;
    real_type const * s;
    real_type const * t;
    s = getVectorPointer(
      arg_in_2, size, CMD "`s` expected to be a real vector"
    );
    t = getVectorPointer(
      arg_in_3, sizet, CMD "`t` expected to be a real vector"
    );
    bool ISO = true;
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string");

    MEX_ASSERT(
      size == sizet || size == 1 || sizet ==1,
      CMD " size(s) = " << size <<
      " must be equal to size(t) = " << sizet <<
      " or size(s|t) == 1"
    );

    mwSize incs = size  == 1 ? 0 : 1;
    mwSize inct = sizet == 1 ? 0 : 1;
    mwSize npts = max(size,sizet);

    #define LOOPXY1 \
      for ( mwSize i = 0; i < npts;  \
            ++i, s += incs, t += inct, pXY += 4 )
    #define LOOPXY2 \
      for ( mwSize i = 0; i < npts; \
            ++i, s += incs, t += inct, ++pTH, ++pK, ++pX, ++pY  )

    if ( nlhs == 1 ) {
      real_type *pXY = createMatrixValue( arg_out_0, 4, size );
      if ( ISO ) {
        LOOPXY1 ptr->evaluate_ISO( *s, *t, pXY[0], pXY[1], pXY[2], pXY[3] );
      } else {
        LOOPXY1 ptr->evaluate_SAE( *s, *t, pXY[0], pXY[1], pXY[2], pXY[3] );
      }
    } else if ( nlhs == 4 ) {
      real_type *pX  = createMatrixValue( arg_out_0, 1, size );
      real_type *pY  = createMatrixValue( arg_out_1, 1, size );
      real_type *pTH = createMatrixValue( arg_out_2, 1, size );
      real_type *pK  = createMatrixValue( arg_out_3, 1, size );
      if ( ISO ) {
        LOOPXY2 ptr->evaluate_ISO( *s, *t, *pTH, *pK, *pX, *pY );
      } else {
        LOOPXY2 ptr->evaluate_SAE( *s, *t, *pTH, *pK, *pX, *pY );
      }
    } else {
      MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
    }
    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD

  } else if ( nrhs == 3 ) {

    #define CMD CMD_BASE "('eval_DDD',OBJ,s): "

    mwSize npts;
    real_type const * s;
    s = getVectorPointer(
      arg_in_2, npts, CMD "`s` expected to be a real vector"
    );

    #define LOOPXY1 \
      for ( mwSize i = 0; i < npts;  \
            ++i, ++s, pXY += 4 )
    #define LOOPXY2 \
      for ( mwSize i = 0; i < npts; \
            ++i, ++s, ++pTH, ++pK, ++pX, ++pY )

    if ( nlhs == 1 ) {
      real_type *pXY = createMatrixValue( arg_out_0, 4, npts );
      LOOPXY1 ptr->evaluate( *s, pXY[0], pXY[1], pXY[2], pXY[3] );
    } else if ( nlhs == 4 ) {
      real_type *pX  = createMatrixValue( arg_out_0, 1, npts );
      real_type *pY  = createMatrixValue( arg_out_1, 1, npts );
      real_type *pTH = createMatrixValue( arg_out_2, 1, npts );
      real_type *pK  = createMatrixValue( arg_out_3, 1, npts );
      LOOPXY2 ptr->evaluate( *s, *pTH, *pK, *pX, *pY);
    } else {
      MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
    }

    #undef LOOPXY1
    #undef LOOPXY2

    #undef CMD
  }
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_theta(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('theta',OBJ,s): "
  MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  mwSize size;
  real_type const * s;
  s = getVectorPointer(
    arg_in_2, size, CMD "`s` expected to be a real vector"
  );

  real_type *theta = createMatrixValue( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i )
    *theta++ = ptr->theta( *s++ );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_theta_D(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('theta_D',OBJ,s): "
  MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  mwSize size;
  real_type const * s;
  s = getVectorPointer(
    arg_in_2, size, CMD "`s` expected to be a real vector"
  );

  real_type *theta = createMatrixValue( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i )
    *theta++ = ptr->theta_D( *s++ );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_theta_DD(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('theta_DD',OBJ,s): "
  MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  mwSize size;
  real_type const * s;
  s = getVectorPointer(
    arg_in_2, size, CMD "`s` expected to be a real vector"
  );

  real_type *theta = createMatrixValue( arg_out_0, 1, size );

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

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('theta_DDD',OBJ,s): "
  MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  mwSize size;
  real_type const * s;
  s = getVectorPointer(
    arg_in_2, size, CMD "`s` expected to be a real vector"
  );

  real_type *theta = createMatrixValue( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i )
    *theta++ = ptr->theta_DDD( *s++ );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_kappa(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('kappa',OBJ,s): "
  MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  mwSize size;
  real_type const * s;
  s = getVectorPointer(
    arg_in_2, size, CMD "`s` expected to be a real vector"
  );

  real_type *kappa = createMatrixValue( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i )
    *kappa++ = ptr->kappa( *s++ );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_kappa_D(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('kappa_D',OBJ,s): "
  MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  mwSize size;
  real_type const * s;
  s = getVectorPointer(
    arg_in_2, size, CMD "`s` expected to be a real vector"
  );

  real_type *kappa_D = createMatrixValue( arg_out_0, 1, size );

  for ( mwSize i = 0; i < size; ++i )
    *kappa_D++ = ptr->kappa_D( *s++ );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_kappa_DD(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('kappa_DD',OBJ,s): "
  MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  mwSize size;
  real_type const * s;
  s = getVectorPointer(
    arg_in_2, size, CMD "`s` expected to be a real vector"
  );

  real_type *kappa_DD = createMatrixValue( arg_out_0, 1, size );

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

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('xyBegin',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = " << nlhs );

  setScalarValue( arg_out_0, ptr->xBegin() );
  setScalarValue( arg_out_1, ptr->yBegin() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_x_begin(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('xBegin',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  setScalarValue( arg_out_0, ptr->xBegin() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_y_begin(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('yBegin',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  setScalarValue( arg_out_0, ptr->yBegin() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_theta_begin(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('thetaBegin',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  setScalarValue( arg_out_0, ptr->thetaBegin() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_kappa_begin(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('kappaBegin',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  setScalarValue( arg_out_0, ptr->kappaBegin() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_xy_end(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('xyEnd',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = " << nlhs );

  setScalarValue( arg_out_0, ptr->xEnd() );
  setScalarValue( arg_out_1, ptr->yEnd() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_x_end(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('xEnd',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  setScalarValue( arg_out_0, ptr->xEnd() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_y_end(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('yEnd',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  setScalarValue( arg_out_0, ptr->yEnd() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_theta_end(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('thetaEnd',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  setScalarValue( arg_out_0, ptr->thetaEnd() );

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

static
void
do_kappa_end(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  G2LIB_CLASS * ptr = convertMat2Ptr<G2LIB_CLASS>(arg_in_1);

  #define CMD CMD_BASE "('kappaEnd',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

  setScalarValue( arg_out_0, ptr->kappaEnd() );

  #undef CMD
}

static
void
do_yesAABBtree(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  #define CMD CMD_BASE "('yesAABBtree',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );

  G2lib::yesAABBtree();

  #undef CMD
}

static
void
do_noAABBtree(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  #define CMD CMD_BASE "('noAABBtree',OBJ): "
  MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
  MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );

  G2lib::noAABBtree();

  #undef CMD
}

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

#define CMD_VIRTUAL_LIST \
CMD_LENGTH,              \
CMD_DELETE,              \
CMD_COPY,                \
CMD_BBOX,                \
CMD_CHANGE_ORIGIN,       \
CMD_TRANSLATE,           \
CMD_ROTATE,              \
CMD_SCALE,               \
CMD_REVERSE,             \
CMD_TRIM,                \
CMD_DISTANCE,            \
CMD_CLOSEST_POINT,       \
CMD_COLLISION,           \
CMD_INTERSECT,           \
CMD_FINDST,              \
CMD_INFO,                \
CMD_EVALUATE,            \
CMD_EVAL,                \
CMD_EVAL_D,              \
CMD_EVAL_DD,             \
CMD_EVAL_DDD,            \
CMD_THETA,               \
CMD_THETA_D,             \
CMD_THETA_DD,            \
CMD_THETA_DDD,           \
CMD_KAPPA,               \
CMD_KAPPA_D,             \
CMD_KAPPA_DD,            \
CMD_XY_BEGIN,            \
CMD_X_BEGIN,             \
CMD_Y_BEGIN,             \
CMD_THETA_BEGIN,         \
CMD_KAPPA_BEGIN,         \
CMD_XY_END,              \
CMD_X_END,               \
CMD_Y_END,               \
CMD_THETA_END,           \
CMD_KAPPA_END,           \
CMD_YES_AABBTREE,        \
CMD_NO_AABBTREE

#define CMD_MAP_LIST                \
{"length",CMD_LENGTH},              \
{"delete",CMD_DELETE},              \
{"copy",CMD_COPY},                  \
{"bbox",CMD_BBOX},                  \
{"changeOrigin",CMD_CHANGE_ORIGIN}, \
{"translate",CMD_TRANSLATE},        \
{"rotate",CMD_ROTATE},              \
{"scale",CMD_SCALE},                \
{"reverse",CMD_REVERSE},            \
{"trim",CMD_TRIM},                  \
{"distance",CMD_DISTANCE},          \
{"closestPoint",CMD_CLOSEST_POINT}, \
{"collision",CMD_COLLISION},        \
{"intersect",CMD_INTERSECT},        \
{"findST",CMD_FINDST},              \
{"info",CMD_INFO},                  \
{"evaluate",CMD_EVALUATE},          \
{"eval",CMD_EVAL},                  \
{"eval_D",CMD_EVAL_D},              \
{"eval_DD",CMD_EVAL_DD},            \
{"eval_DDD",CMD_EVAL_DDD},          \
{"theta",CMD_THETA},                \
{"theta_D",CMD_THETA_D},            \
{"theta_DD",CMD_THETA_DD},          \
{"theta_DDD",CMD_THETA_DDD},        \
{"kappa",CMD_KAPPA},                \
{"kappa_D",CMD_KAPPA_D},            \
{"kappa_DD",CMD_KAPPA_DD},          \
{"xyBegin",CMD_XY_BEGIN},           \
{"xBegin",CMD_X_BEGIN},             \
{"yBegin",CMD_Y_BEGIN},             \
{"thetaBegin",CMD_THETA_BEGIN},     \
{"kappaBegin",CMD_KAPPA_BEGIN},     \
{"xyEnd",CMD_XY_END},               \
{"xEnd",CMD_X_END},                 \
{"yEnd",CMD_Y_END},                 \
{"thetaEnd",CMD_THETA_END},         \
{"kappaEnd",CMD_KAPPA_END},         \
{"yesAABBtree",CMD_YES_AABBTREE},   \
{"noAABBtree",CMD_NO_AABBTREE}

#define CMD_CASE_LIST                         \
case CMD_LENGTH:                              \
  do_length( nlhs, plhs, nrhs, prhs );        \
  break;                                      \
case CMD_DELETE:                              \
  do_delete( nlhs, plhs, nrhs, prhs );        \
  break;                                      \
case CMD_COPY:                                \
  do_copy( nlhs, plhs, nrhs, prhs );          \
  break;                                      \
case CMD_BBOX:                                \
  do_bbox( nlhs, plhs, nrhs, prhs );          \
  break;                                      \
case CMD_CHANGE_ORIGIN:                       \
  do_change_origin( nlhs, plhs, nrhs, prhs ); \
  break;                                      \
case CMD_TRANSLATE:                           \
  do_translate( nlhs, plhs, nrhs, prhs );     \
  break;                                      \
case CMD_ROTATE:                              \
  do_rotate( nlhs, plhs, nrhs, prhs );        \
  break;                                      \
case CMD_SCALE:                               \
  do_scale( nlhs, plhs, nrhs, prhs );         \
  break;                                      \
case CMD_REVERSE:                             \
  do_reverse( nlhs, plhs, nrhs, prhs );       \
  break;                                      \
case CMD_TRIM:                                \
  do_trim( nlhs, plhs, nrhs, prhs );          \
  break;                                      \
case CMD_DISTANCE:                            \
  do_distance( nlhs, plhs, nrhs, prhs );      \
  break;                                      \
case CMD_CLOSEST_POINT:                       \
  do_closestPoint( nlhs, plhs, nrhs, prhs );  \
  break;                                      \
case CMD_COLLISION:                           \
  do_collision( nlhs, plhs, nrhs, prhs );     \
  break;                                      \
case CMD_INTERSECT:                           \
  do_intersect( nlhs, plhs, nrhs, prhs );     \
  break;                                      \
case CMD_FINDST:                              \
  do_findST( nlhs, plhs, nrhs, prhs );        \
  break;                                      \
case CMD_INFO:                                \
  do_info( nlhs, plhs, nrhs, prhs );          \
  break;                                      \
case CMD_EVALUATE:                            \
  do_evaluate( nlhs, plhs, nrhs, prhs );      \
  break;                                      \
case CMD_EVAL:                                \
  do_eval( nlhs, plhs, nrhs, prhs );          \
  break;                                      \
case CMD_EVAL_D:                              \
  do_eval_D( nlhs, plhs, nrhs, prhs );        \
  break;                                      \
case CMD_EVAL_DD:                             \
  do_eval_DD( nlhs, plhs, nrhs, prhs );       \
  break;                                      \
case CMD_EVAL_DDD:                            \
  do_eval_DDD( nlhs, plhs, nrhs, prhs );      \
  break;                                      \
case CMD_THETA:                               \
  do_theta( nlhs, plhs, nrhs, prhs );         \
  break;                                      \
case CMD_THETA_D:                             \
  do_theta_D( nlhs, plhs, nrhs, prhs );       \
  break;                                      \
case CMD_THETA_DD:                            \
  do_theta_DD( nlhs, plhs, nrhs, prhs );      \
  break;                                      \
case CMD_THETA_DDD:                           \
  do_theta_DDD( nlhs, plhs, nrhs, prhs );     \
  break;                                      \
case CMD_KAPPA:                               \
  do_kappa( nlhs, plhs, nrhs, prhs );         \
  break;                                      \
case CMD_KAPPA_D:                             \
  do_kappa_D( nlhs, plhs, nrhs, prhs );       \
  break;                                      \
case CMD_KAPPA_DD:                            \
  do_kappa_DD( nlhs, plhs, nrhs, prhs );      \
  break;                                      \
case CMD_XY_BEGIN:                            \
  do_xy_begin( nlhs, plhs, nrhs, prhs );      \
  break;                                      \
case CMD_X_BEGIN:                             \
  do_x_begin( nlhs, plhs, nrhs, prhs );       \
  break;                                      \
case CMD_Y_BEGIN:                             \
  do_y_begin( nlhs, plhs, nrhs, prhs );       \
  break;                                      \
case CMD_THETA_BEGIN:                         \
  do_theta_begin( nlhs, plhs, nrhs, prhs );   \
  break;                                      \
case CMD_KAPPA_BEGIN:                         \
  do_kappa_begin( nlhs, plhs, nrhs, prhs );   \
  break;                                      \
case CMD_XY_END:                              \
  do_xy_end( nlhs, plhs, nrhs, prhs );        \
  break;                                      \
case CMD_X_END:                               \
  do_x_end( nlhs, plhs, nrhs, prhs );         \
  break;                                      \
case CMD_Y_END:                               \
  do_y_end( nlhs, plhs, nrhs, prhs );         \
  break;                                      \
case CMD_THETA_END:                           \
  do_theta_end( nlhs, plhs, nrhs, prhs );     \
  break;                                      \
case CMD_KAPPA_END:                           \
  do_kappa_end( nlhs, plhs, nrhs, prhs );     \
  break;                                      \
case CMD_YES_AABBTREE:                        \
  do_yesAABBtree( nlhs, plhs, nrhs, prhs );   \
  break;                                      \
case CMD_NO_AABBTREE:                         \
  do_noAABBtree( nlhs, plhs, nrhs, prhs );    \
  break
