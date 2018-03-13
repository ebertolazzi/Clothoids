/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Circle.hh"
#include "mex_utils.hh"

#include <vector>

#define MEX_ERROR_MESSAGE \
"==========================================================================\n" \
"Compute cicle arc\n" \
"\n" \
"USAGE:\n" \
"\n" \
"  OBJ = CircleMexWrapper( 'new', x0, y0, theta0, kur, L ) ;\n" \
"  OBJ = CircleMexWrapper( 'new', p0, p1, p2 ) ;\n" \
"  OBJ = CircleMexWrapper( 'new', p0, theta0, p1 ) ;\n" \
"\n" \
"  CircleMexWrapper( 'delete', OBJ ) ;\n" \
"\n" \
"  CircleMexWrapper( 'build', OBJ, x0, y0, theta0, kur, L ) ;\n" \
"  CircleMexWrapper( 'build', OBJ, p0, p1, p2 ) ;\n" \
"  CircleMexWrapper( 'build', OBJ, p0, theta0, p1 ) ;\n" \
"\n" \
"  CircleMexWrapper( 'changeOrigin', OBJ, x0, y0 ) ;\n" \
"  CircleMexWrapper( 'translate', OBJ, tx, ty ) ;\n" \
"  CircleMexWrapper( 'trim', OBJ, smin, smax ) ;\n" \
"  CircleMexWrapper( 'rotate', OBJ, angle, cx, cy ) ;\n" \
"  CircleMexWrapper( 'scale', OBJ, scale ) ;\n" \
"  CircleMexWrapper( 'reverse', OBJ ) ;\n" \
"\n" \
"  [d,s] = CircleMexWrapper( 'distance', OBJ, x, y ) ;\n" \
"  [d,s] = CircleMexWrapper( 'distance', OBJ, p ) ;\n" \
"\n" \
"  [p0,p1,p2,ok] = CircleMexWrapper( 'bbTriangle', OBJ ) ;\n" \
"\n" \
"  burbs = CircleMexWrapper( 'to_nurbs', OBJ ) ;\n" \
"\n" \
"  res = CircleMexWrapper( 'getX0', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'getY0', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'getTheta0', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'getKappa', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'getSmin', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'getSmax', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'length', OBJ ) ;\n" \
"\n" \
"  [X,Y] = CircleMexWrapper( 'eval', OBJ, s ) ;\n" \
"  [X,Y] = CircleMexWrapper( 'eval_D', OBJ, s ) ;\n" \
"  [X,Y] = CircleMexWrapper( 'eval_DD', OBJ, s ) ;\n" \
"  [X,Y] = CircleMexWrapper( 'eval_DDD', OBJ, s ) ;\n" \
"\n" \
"  [X,Y] = CircleMexWrapper( 'theta', OBJ, s ) ;\n" \
"\n" \
"  nurbs = CircleMexWrapper( 'to_nurbs', OBJ ) ;\n" \
"\n" \
"%==========================================================================%\n" \
"%                                                                          %\n" \
"%  Autor: Enrico Bertolazzi                                                %\n" \
"%         Department of Industrial Engineering                             %\n" \
"%         University of Trento                                             %\n" \
"%         enrico.bertolazzi@unitn.it                                       %\n" \
"%                                                                          %\n" \
"%==========================================================================%\n"

namespace G2lib {

  using namespace std;

  static
  CircleArc *
  DATA_NEW( mxArray * & mx_id ) {
    CircleArc * ptr = new CircleArc();
    mx_id = convertPtr2Mat<CircleArc>(ptr);
    return ptr ;
  }

  static
  inline
  CircleArc *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<CircleArc>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<CircleArc>(mx_id);
  }

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    // the first argument must be a string
    if ( nrhs == 0 ) {
      mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    }

    try {

      MEX_ASSERT( mxIsChar(arg_in_0), "CircleArc(...): First argument must be a string" ) ;
      string cmd = mxArrayToString(arg_in_0) ;
      mwSize size0, size1, size2 ;

      bool do_new = cmd == "new" ;
      CircleArc * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new || cmd == "build" ) {

        indexType kk = do_new ? 0 : 1 ;

        if ( do_new ) {
          MEX_ASSERT( nlhs == 1, "CircleArc, expected 1 output" );
        } else {
          MEX_ASSERT( nlhs == 0, "CircleArc, expected no output" );
        }

        if ( nrhs == 6+kk ) {

          valueType x0     = getScalarValue( prhs[1+kk], "CircleArc(x0,y0,theta0,k0,...): `x0` expected to be a real scalar" );
          valueType y0     = getScalarValue( prhs[2+kk], "CircleArc(x0,y0,theta0,k0,...): `y0` expected to be a real scalar" );
          valueType theta0 = getScalarValue( prhs[3+kk], "CircleArc(x0,y0,theta0,k0,...): `theta0` expected to be a real scalar" );
          valueType k0     = getScalarValue( prhs[4+kk], "CircleArc(x0,y0,theta0,k0,...): `k0` expected to be a real scalar" );
          valueType L      = getScalarValue( prhs[5+kk], "CircleArc(x0,y0,theta0,k0,L): `L` expected to be a real scalar" );

          ptr->build( x0, y0, theta0, k0, L );

        } else if ( nrhs == 4+kk ) {

          valueType const * p0 = getVectorPointer( prhs[1+kk], size0, "CircleArc(p0,p1 or theta0,p2): `p0` expected to be a real vector" );
          valueType const * p1 = getVectorPointer( prhs[2+kk], size1, "CircleArc(p0,p1 or theta0,p2): `p1` expected to be a real vector" );
          valueType const * p2 = getVectorPointer( prhs[3+kk], size2, "CircleArc(p0,p1 or theta0,p2): `p2` expected to be a real vector" );

          MEX_ASSERT( size0 == 2 && size2 == 2 && size1 > 0 && size1 <= 2,
                      "CircleArc: bad dimension size(p0) = " << size0 <<
                      ", size(p1) = " << size1 << ", size(p2) = " << size2 ) ;

          if ( size1 == 1 ) ptr->build_G1( p0[0], p0[1], p1[0] /* theta0 */, p2[0], p2[1] ) ;
          else              ptr->build_3P( p0[0], p0[1], p1[0], p1[1],       p2[0], p2[1] ) ;

        } else if ( nrhs == 1 ) {
          // nothing to do
        } else {
          MEX_ASSERT(false, "CircleArc, expected " << 4+kk << ", " << 7+kk << " or " << 8+kk << " inputs, nrhs = " << nrhs );
        }

        plhs[0] = convertPtr2Mat<CircleArc>(ptr);

      } else if ( cmd == "delete" ) {

        MEX_ASSERT(nrhs == 2, "CircleArc('delete',OBJ): expected 2 inputs");
        MEX_ASSERT(nlhs == 0, "CircleArc('delete',OBJ): expected no output");
        // Destroy the C++ object
        DATA_DELETE( arg_in_1 ) ;

      } else if ( cmd == "changeOrigin" ) {

        MEX_ASSERT(nrhs == 4, "CircleArc('changeOrigin',OBJ,x0,y0): expected 4 inputs");
        MEX_ASSERT(nlhs == 0, "CircleArc('changeOrigin',OBJ,x0,y0): expected no output");

        valueType new_x0 = getScalarValue( arg_in_2, "CircleArc('changeOrigin',OBJ,x0,y0): `x0` expected to be a real scalar" );
        valueType new_y0 = getScalarValue( arg_in_3, "CircleArc('changeOrigin',OBJ,x0,y0): `y0` expected to be a real scalar" );

        ptr->changeOrigin( new_x0, new_y0 );

      } else if ( cmd == "translate" ) {

        MEX_ASSERT(nrhs == 4, "CircleArc('translate',OBJ,tx,ty): expected 4 inputs");
        MEX_ASSERT(nlhs == 0, "CircleArc('translate',OBJ,tx,ty): expected no output");

        valueType tx = getScalarValue( arg_in_2, "CircleArc('translate',OBJ,tx,ty): `tx` expected to be a real scalar" );
        valueType ty = getScalarValue( arg_in_3, "CircleArc('translate',OBJ,tx,ty): `ty` expected to be a real scalar" );

        ptr->translate( tx, ty );

      } else if ( cmd == "changeCurvilinearOrigin" ) {

        MEX_ASSERT(nrhs == 4, "CircleArc('changeCurvilinearOrigin',OBJ,s0,L): expected 4 inputs");

        valueType s0 = getScalarValue(arg_in_2,"CircleArc('changeCurvilinearOrigin',OBJ,s0,L): Error in reading s0") ;
        valueType L  = getScalarValue(arg_in_3,"CircleArc('changeCurvilinearOrigin',OBJ,s0,L): Error in reading L") ;
        ptr->changeCurvilinearOrigin(s0,L);

      } else if ( cmd == "rotate" ) {

        MEX_ASSERT(nrhs == 5, "CircleArc('rotate',OBJ,angle,cx,cy): expected 5 inputs");
        MEX_ASSERT(nlhs == 0, "CircleArc('rotate',OBJ,angle,cx,cy): expected no output");

        valueType angle = getScalarValue( arg_in_2, "CircleArc('rotate',OBJ,angle,cx,cy): `angle` expected to be a real scalar" );
        valueType cx    = getScalarValue( arg_in_3, "CircleArc('rotate',OBJ,angle,cx,cy): `cx` expected to be a real scalar" );
        valueType cy    = getScalarValue( arg_in_4, "CircleArc('rotate',OBJ,angle,cx,cy): `cy` expected to be a real scalar" );

        ptr->rotate( angle, cx, cy );

      } else if ( cmd == "scale" ) {

        MEX_ASSERT(nrhs == 3, "CircleArc('scale',OBJ,scale): expected 3 inputs");
        MEX_ASSERT(nlhs == 0, "CircleArc('scale',OBJ,scale): expected no output");

        valueType sc = getScalarValue( arg_in_2, "CircleArc('scale',OBJ,scale): `scale` expected to be a real scalar" );
        ptr->scale( sc );

      } else if ( cmd == "reverse" ) {

        MEX_ASSERT(nrhs == 2, "CircleArc('reverse',OBJ): expected 2 inputs");
        MEX_ASSERT(nlhs == 0, "CircleArc('reverse',OBJ):expected no output");
        ptr->reverse();

      } else if ( cmd == "trim" ) {

        MEX_ASSERT(nrhs == 4, "CircleArc('trim',OBJ,s_begin,s_end): expected 4 inputs");
        MEX_ASSERT(nlhs == 0, "CircleArc('trim',OBJ,s_begin,s_end): expected no output");

        valueType s_begin = getScalarValue( arg_in_2, "CircleArc('trim',OBJ,s_begin,s_end): `s_begin` expected to be a real scalar" );
        valueType s_end   = getScalarValue( arg_in_3, "CircleArc('trim',OBJ,s_begin,s_end): `s_end` expected to be a real scalar" );

        ptr->trim( s_begin, s_end );

      } else if ( cmd == "theta" ) {

        MEX_ASSERT(nrhs == 3, "CircleArc('theta',OBJ,s): expected 3 inputs");
        MEX_ASSERT(nlhs == 1, "CircleArc('theta',OBJ,s): expected 1 output");

        valueType s = getScalarValue( arg_in_2, "CircleArc('theta',OBJ,s): `s` expected to be a real scalar" );
        setScalarValue( arg_out_0, ptr->theta( s ) ) ;

      } else if ( cmd == "distance" ) {
        valueType pnt[2] ;
        valueType const * p = pnt ;

        switch ( nrhs ) {
        case 4:
          pnt[0] = getScalarValue( arg_in_2, "CircleArc('distance',OBJ,x,y): `x` expected to be a real scalar" );
          pnt[1] = getScalarValue( arg_in_3, "CircleArc('distance',OBJ,x,y): `y` expected to be a real scalar" );
          break ;
        case 3:
          p = getVectorPointer( arg_in_2, size0, "CircleArc('distance',OBJ,pnt): `pnt` expected to be a real vector" );
          break ;

        default:
          MEX_ASSERT( false, "CircleArc('distance',OBJ,x,y): expected 3 or 4 inputs");
          break ;
        }

        valueType ss ;
        valueType dst = ptr->distance( p[0], p[1], ss );

        MEX_ASSERT(nlhs >= 1 && nlhs <= 2, "CircleArc('distance',OBJ,pnt): expected 1,2 or 3 output");

        if ( nlhs > 0 ) setScalarValue( arg_out_0, dst ) ;
        if ( nlhs > 1 ) setScalarValue( arg_out_1, ss ) ;

      } else if ( cmd == "bbTriangle" ) {

        MEX_ASSERT(nrhs == 2, "CircleArc('bbTriangle',OBJ): expected 2 inputs");
        MEX_ASSERT(nlhs == 4, "CircleArc('bbTriangle',OBJ): expected 4 output");

        double * p0 = createMatrixValue( arg_out_0, 1, 2 );
        double * p1 = createMatrixValue( arg_out_1, 1, 2 );
        double * p2 = createMatrixValue( arg_out_2, 1, 2 );

        bool ok = ptr->bbTriangle( p0, p1, p2 );

        arg_out_3 = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
        *static_cast<mxLogical*>(mxGetData(arg_out_3)) = ok ;

      } else if ( cmd == "to_nurbs" ) {

        valueType knots[12], Poly[9][3] ;
        indexType npts = ptr->toNURBS( knots, Poly ); // npt + 2

        static char const * fieldnames[] = { "form", "order", "dim", "number", "knots", "coefs" } ;
        arg_out_0 = mxCreateStructMatrix(1,1,6,fieldnames);
        mxArray * mx_knots = mxCreateDoubleMatrix(1,npts+2,mxREAL);
        mxArray * mx_Poly  = mxCreateDoubleMatrix(3,npts,mxREAL);

        mxSetFieldByNumber( arg_out_0, 0, 0, mxCreateString("rB") );
        mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(3) );
        mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar(2) );
        mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar(npts) );
        mxSetFieldByNumber( arg_out_0, 0, 4, mx_knots );
        mxSetFieldByNumber( arg_out_0, 0, 5, mx_Poly );

        double *kb = mxGetPr(mx_knots) ;
        for ( indexType i = 0 ; i < npts+2 ; ++i ) *kb++ = knots[i] ;

        double *pr = mxGetPr(mx_Poly) ;
        for ( indexType i = 0 ; i < npts ; ++i ) {
          *pr++ = Poly[i][0] ;
          *pr++ = Poly[i][1] ;
          *pr++ = Poly[i][2] ;
        }
      } else {
        if ( nrhs == 3 ) {
          mwSize size ;
          double const * s = getVectorPointer( arg_in_2, size, "CircleArc('eval*',s): `s` expected to be a real vector" ) ;
          double *pX = createMatrixValue( arg_out_0, 1,size );
          double *pY = createMatrixValue( arg_out_1, 1,size );
          if ( cmd == "eval" ) {
            for ( mwSize i = 0 ; i < size ; ++i, ++s, ++pX, ++pY )
              ptr->eval( *s, *pX, *pY ) ;
          } else if ( cmd == "eval_D" ) {
            for ( mwSize i = 0 ; i < size ; ++i, ++s, ++pX, ++pY )
              ptr->eval_D( *s, *pX, *pY ) ;
          } else if ( cmd == "eval_DD" ) {
            for ( mwSize i = 0 ; i < size ; ++i, ++s, ++pX, ++pY )
              ptr->eval_DD( *s, *pX, *pY ) ;
          } else if ( cmd == "eval_DDD" ) {
            for ( mwSize i = 0 ; i < size ; ++i, ++s, ++pX, ++pY )
              ptr->eval_DDD( *s, *pX, *pY ) ;
          } else {
            MEX_ASSERT(false, "Unknown command: " << cmd );
          }
        } else if ( nrhs == 2 ) {
          if      ( cmd == "getX0"     ) setScalarValue( arg_out_0, ptr->getX0());
          else if ( cmd == "getY0"     ) setScalarValue( arg_out_0, ptr->getY0());
          else if ( cmd == "getTheta0" ) setScalarValue( arg_out_0, ptr->getTheta0());
          else if ( cmd == "getKappa"  ) setScalarValue( arg_out_0, ptr->getKappa());
          else if ( cmd == "length"    ) setScalarValue( arg_out_0, ptr->totalLength());
          else {
            MEX_ASSERT(false, "Unknown command: " << cmd );
          }
        } else {
          MEX_ASSERT(false, "Unknown command: " << cmd );
        }
      }

    } catch ( exception const & e ) {
    	mexErrMsgTxt(e.what()) ;
    } catch (...) {
    	mexErrMsgTxt("circleArc failed\n") ;
    }
  }
}
