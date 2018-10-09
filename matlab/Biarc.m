classdef Biarc < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    function self = Biarc( varargin )
      %% Create a new C++ class instance for the clothoid arc object
      % Usage:
      %    ref = Biarc()
      %    ref = Biarc( x0, y0, theta0, x1, y1, theta1 )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the clothoid at initial point
      %    x1, y1: coordinate of final point
      %    theta1: orientation of the clothoid at final point
      %
      %  On output:
      %    ref: reference handle to the object instance
      self.objectHandle = BiarcMexWrapper( 'new' );
      if nargin > 0
        BiarcMexWrapper( 'build', self.objectHandle, varargin{:} );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      %% Destroy the C++ class instance
      BiarcMexWrapper( 'delete', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = build( self, x0, y0, theta0, x1, y1, theta1 )
      % Build the interpolating G1 biarc
      %
      % Usage:
      %    ref.build_G1( x0, y0, theta0, x1, y1, theta1 )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the clothoid at initial point
      %    x1, y1: coordinate of final point
      %    theta1: orientation of the clothoid at final point
      %
      ok = BiarcMexWrapper( 'build', ...
                            self.objectHandle, ...
                            x0, y0, theta0, ...
                            x1, y1, theta1 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function str = is_type( ~ )
      str = 'BiArc' ;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = build_3P( self, varargin )
      % Build the interpolating biarc by 3 points
      %
      % Usage:
      %    ref.build_3P( x0, y0, x1, y1, x2, y2 )
      %    ref.build_3P( [x0, y0], [x1, y1], [x2, y2] )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    x1, y1: coordinate of middle point
      %    x2, y2: coordinate of final point
      % alternative
      %    p0: coordinate of initial point
      %    p1: coordinate of middle point
      %    p2: coordinate of final point
      %
      ok = BiarcMexWrapper( 'build_3P', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = evaluate( self, s )
      % evaluate the curve at curvilinear abscissa `s`
      %
      % Usage:
      %    [x,y] = ref.eval( s )
      %    [x,y,theta,kappa] = ref.eval( s )
      %
      % On input:
      %    s: curvilinear coordinates where to evaluate the curve
      %       (scalar or vector)
      %
      % On output:
      %    x, y:  coordinates of the curve
      %    theta: orientation of the curve
      %    kappa: curvature of the curve
      %
      [ varargout{1:nargout} ] = ...
        BiarcMexWrapper( 'evaluate', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [x,y] = eval( self, varargin )
      [ x, y ] = BiarcMexWrapper( 'eval', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [x_D,y_D] = eval_D( self, varargin )
      [ x_D, y_D ] = ...
        BiarcMexWrapper( 'eval_D', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [x_DD,y_DD] = eval_DD( self, varargin )
      [ x_DD, y_DD ] = ...
        BiarcMexWrapper( 'eval_DD', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [x_DDD,y_DDD] = eval_DDD( self, varargin )
      [ x_DDD, y_DDD ] = ...
        BiarcMexWrapper( 'eval_DDD', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [X,Y,S,DST] = closestPoint( self, qx, qy )
      [ X, Y, S, DST ] = ...
        BiarcMexWrapper( 'closestPoint', self.objectHandle, qx, qy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [DST,S] = distance( self, varargin )
      % eval the angle of the circle curve at curvilinear abscissa `s`
      [ DST, S ] = ...
        BiarcMexWrapper( 'distance', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [X,Y,S,DST] = closestPointBySample( self, qx, qy, ds )
      [ X, Y, S, DST ] = ...
         BiarcMexWrapper( 'closestPointBySample', ...
                          self.objectHandle, qx, qy, ds );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [DST,S] = distanceBySample( self, qx, qy, ds )
      [ DST, S ] = ...
        BiarcMexWrapper( 'distanceBySample', ...
                         self.objectHandle, qx, qy, ds );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = xBegin0( self )
      res = BiarcMexWrapper( 'xBegin0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = xEnd0( self )
      res = BiarcMexWrapper( 'xEnd0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = xBegin1( self )
      res = BiarcMexWrapper( 'xBegin1', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = xEnd1( self )
      res = BiarcMexWrapper( 'xEnd1', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = yBegin0( self )
      res = BiarcMexWrapper( 'yBegin0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = yEnd0( self )
      res = BiarcMexWrapper( 'yEnd0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = yBegin1( self )
      res = BiarcMexWrapper( 'yBegin1', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = yEnd1( self )
      res = BiarcMexWrapper( 'yEnd1', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = thetaBegin0( self )
      res = BiarcMexWrapper( 'thetaBegin0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = thetaBegin1( self )
      res = BiarcMexWrapper( 'thetaBegin1', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = thetaEnd0( self )
      res = BiarcMexWrapper( 'thetaEnd0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = thetaEnd1( self )
      res = BiarcMexWrapper( 'thetaEnd1', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = kappa0( self )
      res = BiarcMexWrapper( 'kappa0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = kappa1( self )
      res = BiarcMexWrapper( 'kappa1', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = length0( self )
      res = BiarcMexWrapper( 'length0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = length1( self )
      res = BiarcMexWrapper( 'length1', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = length( self )
      res = BiarcMexWrapper( 'length', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function rotate( self, angle, cx, cy )
      % rotate the clothoid curve by angle respect to the centre `(cx,cy)`
      %
      % Usage:
      %    ref.rotate(angle, cx, cy)
      %
      % On input:
      %    angle: the angle of rotation
      %    cx, cy: coordinates of the centre of rotation
      %
      BiarcMexWrapper( 'rotate', self.objectHandle, angle, cx, cy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function translate( self, tx, ty )
      % translate the clothoid curve by `(tx,ty)`
      %
      % Usage:
      %    ref.translate(tx, ty)
      %
      % On input:
      %    tx, ty: horizontal and vertical translation
      %
      BiarcMexWrapper( 'translate', self.objectHandle, tx, ty );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function changeOrigin( self, newX0, newY0 )
      % move the origin of the clothoid to `(newX0, newY0)`
      %
      % Usage:
      %    ref.changeOrigin(newX0, newY0)
      %
      % On input:
      %    newX0, newY0: new coordinates of initial point
      %
      BiarcMexWrapper( 'changeOrigin', self.objectHandle, newX0, newY0 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function scale( self, s )
      % scale clothoid by `sc` factor
      %
      % Usage:
      %    ref.scale(newX0, newY0)
      %
      % On input:
      %    newX0, newY0: new coordinates of initial point

      BiarcMexWrapper( 'scale', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function reverse( self )
      % reverse the orientation of the clothoid curve
      % Usage:
      %    ref.reverse()
      %
      BiarcMexWrapper( 'reverse', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ C0, C1 ] = getCircles( self )
      C0 = CircleArc() ;
      C1 = CircleArc() ;
      C0.build( self.xBegin0(), self.yBegin0(), ...
                self.thetaBegin0(), self.kappa0(), self.length0());
      C1.build( self.xBegin1(), self.yBegin1(), ...
                self.thetaBegin1(), self.kappa1(), self.length1());
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ arc0, arc1 ] = to_nurbs( self )
      % Usage:
      %    ref.to_nurbs()
      %
      [ arc0, arc1 ] = BiarcMexWrapper( 'to_nurbs', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ s, t ] = find_coord( self, x, y )
      [ s, t ] = BiarcMexWrapper( 'findST', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function info( self )
      BiarcMexWrapper( 'info', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, npts, varargin )
      if nargin<2
        npts = 64 ;
      end
      if nargin>2
        fmt1 = varargin{1};
      else
        fmt1 = {'Color','blue','Linewidth',2} ;
      end
      if nargin>3
        fmt2 = varargin{2};
      else
        fmt2 = {'Color','red','Linewidth',2} ;
      end
      [C0,C1] = self.getCircles();
      C0.plot(npts,fmt1) ;
      hold on ;
      C1.plot(npts,fmt2) ;
    end
  end
end
