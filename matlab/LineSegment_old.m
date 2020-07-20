classdef LineSegment < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    function self = LineSegment( varargin )
      %% Create a new C++ class instance for the Segment object
      %
      % Usage:
      %    (1) ref = LineSegment()
      %    (2) ref = LineSegment( x0, y0, theta0, L )
      %    (3) ref = LineSegment( x0, y0, theta0, smin, smax )
      %    (4) ref = LineSegment( p0, p1 )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the circle at initial point
      %    L:      length of curve from initial to final point
      %    smin:   initial curvilinear coordinate of the curve
      %    smax:   final curvilinear coordinate of the curve
      %    p0:     2D point
      %    p1:     2D point
      %
      %    (1) empty segment
      %    (2) line segment passing from (x0,y0) at angle theta0
      %    (3) line segment as in (2) with intial and final curvilinear coordinate respect to (x0,y0)
      %    (4) segment passing from 2 points
      %  On output:
      %    ref: reference handle to the object instance
      %
      self.objectHandle = LineSegmentMexWrapper('new', varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
      %% Destroy the C++ class instance
      LineSegmentMexWrapper('delete', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function obj = obj_handle( self )
      obj = self.objectHandle;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function str = is_type( ~ )
      str = 'LineSegment';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function copy( C )
      LineSegmentMexWrapper( 'copy', self.objectHandle, C.obj_handle() );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function build( self, varargin )
      %
      % Usage:
      %    (1) ref.build( x0, y0, theta0, k0, L )
      %    (3) ref.build( p0, p1 )
      %
      LineSegmentMexWrapper('build', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function translate( self, tx, ty )
      % move the object by `(tx,ty)`
      LineSegmentMexWrapper('translate', self.objectHandle, tx, ty );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function trim( self, smin, smax )
      % trim circle curve to the corresponding curvilinear parameters
      LineSegmentMexWrapper('trim', self.objectHandle, smin, smax );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function rotate( self, angle, cx, cy )
      % rotate the circle by angle with center of rotation `(cx,cy)`
      % Usage:
      %    ref.rotate(angle, cx, cy)
      %
      % On input:
      %    angle: the angle of rotation
      %    cx, cy: coordinates of the centre of rotation
      %
      LineSegmentMexWrapper( 'rotate', self.objectHandle, angle, cx, cy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function changeOrigin( self, newX0, newY0 )
      % move the origin of the circle to `(newX0, newY0)`
      % Usage:
      %    ref.changeOrigin(newX0, newY0)
      %
      % On input:
      %    newX0, newY0: new coordinates of initial point
      %
      LineSegmentMexWrapper( 'changeOrigin', ...
                             self.objectHandle, newX0, newY0 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function nurbs = to_nurbs( self )
      % return a nurbs representation of the circle arc
      nurbs = LineSegmentMexWrapper( 'to_nurbs', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval( self, varargin )
      % eval the circle at curvilinear abscissa `s`
      % Usage:
      %    [x,y] = ref.eval( s )
      %
      % On input:
      %    s: curvilinear coordinates where to evaluate the curve
      %       (scalar or vector)
      %
      % On output:
      %    x, y:  coordinates of the curve
      %
      [ varargout{1:nargout} ] = ...
        LineSegmentMexWrapper( 'eval', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_D( self, varargin )
      % eval the circle derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        LineSegmentMexWrapper( 'eval_D', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DD( self, varargin )
      % eval the circle second derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        LineSegmentMexWrapper( 'eval_DD', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DDD( self, varargin )
      % eval the circle third derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        LineSegmentMexWrapper( 'eval_DDD', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function X0 = xBegin( self )
      X0 = LineSegmentMexWrapper('xBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Y0 = yBegin( self )
      Y0 = LineSegmentMexWrapper('yBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th0 = thetaBegin( self )
      th0 = LineSegmentMexWrapper('thetaBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function X1 = xEnd( self )
      X1 = LineSegmentMexWrapper('xEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Y1 = yEnd( self )
      Y1 = LineSegmentMexWrapper('yEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th1 = thetaEnd( self )
      th1 = LineSegmentMexWrapper('thetaEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = length( self, varargin )
      res = LineSegmentMexWrapper('length', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [p1,p2] = points( self )
      [p1,p2] = LineSegmentMexWrapper('points', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [d,s] = distance( self, x, y, varargin )
      [d,s] = LineSegmentMexWrapper( 'distance', self.objectHandle, ...
                                     x, y, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [s,idx] = projection( self, x, y, varargin )
      [d,s] = LineSegmentMexWrapper( 'projection', self.objectHandle, ...
                                     x, y, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = collision( self, OBJ, varargin )
      [d,s] = LineSegmentMexWrapper( 'collision', ...
                                     self.objectHandle, ...
                                     OBJ.obj_handle(), OBJ.is_type(), ...
                                     varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [s1,s2] = intersect( self, OBJ, varargin )
      [s1,s2] = LineSegmentMexWrapper( 'intersect', ...
                                        self.objectHandle, ...
                                        OBJ.obj_handle(), OBJ.is_type(), ...
                                        varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function info( self )
      LineSegmentMexWrapper( 'info', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [s,t] = find_coord( self, x, y )
      [s,t] = LineSegmentMexWrapper( 'findST', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, varargin )
      [ p1, p2 ] = self.points();
      plot( [ p1(1), p2(1) ], [ p1(2), p2(2) ], varargin{:} );
    end
  end
end
