classdef CircleArc < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    function self = CircleArc( varargin )
      %% Create a new C++ class instance for the circle arc object
      %
      % Usage:
      %    (1) ref = CircleArc()
      %    (2) ref = CircleArc( x0, y0, theta0, k0, L )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the circle at initial point
      %    k0:     curvature of the circle at initial point
      %    L:      length of curve from initial to final point
      %
      %    (1) empty circle
      %    (2) circle passing from (x0,y0) at angle theta0 with curvature and length
      %  On output:
      %    ref: reference handle to the object instance
      %
      self.objectHandle = CircleArcMexWrapper( 'new', varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
      %% Destroy the C++ class instance
      CircleArcMexWrapper( 'delete', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function obj = obj_handle( self )
      obj = self.objectHandle ;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function str = is_type( ~ )
      str = 'CircleArc' ;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function copy( C )
      CircleArcMexWrapper( 'copy', self.objectHandle, C.obj_handle() );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function build( self, x0, y0, theta0, k0, L )
      %
      % Build the circle from known parameters
      %
      % Usage:
      %    (1) ref.build( x0, y0, theta0, k0, L )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the circle at initial point
      %    k0:     curvature of the circle at initial point
      %    L:      length of curve from initial to final point
      %
      %    (1) circle passing from (x0,y0) at angle theta0 with curvature and length
      CircleArcMexWrapper( 'build', ...
                           self.objectHandle, x0, y0, theta0, k0, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = build_G1( self, varargin )
      %
      % Build the circle from known parameters
      %
      % Usage:
      %    (1) ref.build( x0, y0, theta0, x1, y1 )
      %    (2) ref.build( p0, theta0, p1 )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the circle at initial point
      %    k0:     curvature of the circle at initial point
      %    L:      length of curve from initial to final point
      %    p0:     2D point
      %    p1:     2D point
      %
      %    (1) circle passing to [x0,y0] and [x1,y1] with angle theta0 at [x0,y0]
      %    (2) circle passing to p0 and p1 with angle theta0 at p0
      %
      ok = CircleArcMexWrapper( 'build_G1', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = build_3P( self, varargin )
      %
      % Build the circle from known parameters
      %
      % Usage:
      %    (1) ref.build( x0, y0, x1, y1, x2, y2 )
      %    (2) ref.build( p0, p1, p2 )
      %
      ok = CircleArcMexWrapper( 'build_3P', ...
                                self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function translate( self, tx, ty )
      % move the object by `(tx,ty)`
      CircleArcMexWrapper( 'translate', self.objectHandle, tx, ty );
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
      CircleArcMexWrapper( 'rotate', self.objectHandle, angle, cx, cy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function scale( self, sc )
      % scale circle by `sc` factor
      CircleArcMexWrapper( 'scale', self.objectHandle, sc );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function trim( self, smin, smax )
      % trim circle curve to the corresponding curvilinear parameters
      CircleArcMexWrapper( 'trim', self.objectHandle, smin, smax );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function changeCurvilinearOrigin( self, s0, L )
      % change the origin of the circle curve to `s0`
      CircleArcMexWrapper( 'changeCurvilinearOrigin', self.objectHandle, s0, L );
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
      CircleArcMexWrapper( 'changeOrigin', self.objectHandle, newX0, newY0 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [xmin,ymin,xmax,ymax] = bbox( self, varargin )
      % return the bounding box triangle of the circle arc
      [xmin,ymin,xmax,ymax] = CircleArcMexWrapper( 'bbox', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [p0,p1,p2,ok] = bbTriangle( self, varargin )
      % return the bounding box triangle of the circle arc
      [p0,p1,p2,ok] = CircleArcMexWrapper( 'bbTriangle', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function nurbs = to_nurbs( self )
      % return a nurbs representation of the circle arc
      nurbs = CircleArcMexWrapper( 'to_nurbs', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval( self, varargin )
      % eval the circle at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        CircleArcMexWrapper( 'eval', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_D( self, varargin )
      % eval the circle derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        CircleArcMexWrapper( 'eval_D', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DD( self, varargin )
      % eval the circle second derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        CircleArcMexWrapper( 'eval_DD', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DDD( self, varargin )
      % eval the circle third derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        CircleArcMexWrapper( 'eval_DDD', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = theta(self, s)
      % eval the angle of the circle curve at curvilinear abscissa `s`
      th = CircleArcMexWrapper( 'theta', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = theta_D(self, s)
      % eval the angle of the circle curve at curvilinear abscissa `s`
      th = CircleArcMexWrapper( 'theta_D', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = theta_DD(self, s)
      % eval the angle of the circle curve at curvilinear abscissa `s`
      th = CircleArcMexWrapper( 'theta_DD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = theta_DDD(self, s)
      % eval the angle of the circle curve at curvilinear abscissa `s`
      th = CircleArcMexWrapper( 'theta_DDD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function X0 = xBegin(self)
      X0 = CircleArcMexWrapper( 'xBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function X1 = xEnd(self)
      X1 = CircleArcMexWrapper( 'xEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Y0 = yBegin(self)
      Y0 = CircleArcMexWrapper( 'yBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Y1 = yEnd(self)
      Y1 = CircleArcMexWrapper( 'yEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th0 = thetaBegin(self)
      th0 = CircleArcMexWrapper( 'thetaBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th1 = thetaEnd(self)
      th1 = CircleArcMexWrapper( 'thetaEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function k = kappa(self)
      k = CircleArcMexWrapper( 'theta_D', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = length( self, varargin )
      res = CircleArcMexWrapper( 'length', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [d,s] = distance( self, x, y, varargin )
      [d,s] = CircleArcMexWrapper( 'distance', self.objectHandle, ...
                                    x, y, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [s,idx] = projection( self, x, y, varargin )
      [d,s] = CircleArcMexWrapper( 'projection', self.objectHandle, ...
                                   x, y, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = collision( self, OBJ, varargin )
      [d,s] = CircleArcMexWrapper( 'collision', ...
                                   self.objectHandle, ...
                                   OBJ.obj_handle(), OBJ.is_type(), ...
                                   varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [s1,s2] = intersect( self, OBJ, varargin )
      [s1,s2] = CircleArcMexWrapper( 'intersect', ...
                                     self.objectHandle, ...
                                     OBJ.obj_handle(), OBJ.is_type(), ...
                                     varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ s, t ] = find_coord( self, x, y )
      [ s, t ] = CircleArcMexWrapper( 'findST', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function info( self )
      CircleArcMexWrapper( 'info', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, npts, fmt )
      if nargin<2
        npts = 128 ;
      end
      if nargin<3
        fmt = {'Color','blue','Linewidth',2};
      end
      L = self.length();
      S = 0:L/npts:L ;
      [X,Y] = self.eval(S);
      plot(X,Y,fmt{:});
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotPolygon( self, varargin )
      arc = self.to_nurbs() ;
      xx  = arc.coefs(1,:)./arc.coefs(3,:);
      yy  = arc.coefs(2,:)./arc.coefs(3,:);
      if nargin > 1
        plot(xx,yy,varargin{:});
      else
        plot(xx,yy,':ok');
      end
    end
  end
end
