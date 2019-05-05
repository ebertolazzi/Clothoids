classdef Biarc < CurveBase
  %% MATLAB class wrapper for the underlying C++ class

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
      self@CurveBase( 'BiarcMexWrapper' );
      self.objectHandle = BiarcMexWrapper( 'new' );
      if nargin > 0
        ok = BiarcMexWrapper( 'build', self.objectHandle, varargin{:} );
        if ~ok
          error('Biarc constructor failed');
        end
      end
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
      str = 'BiArc';
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
    function res = xMiddle( self )
      res = BiarcMexWrapper( 'xMiddle', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = yMiddle( self )
      res = BiarcMexWrapper( 'yMiddle', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = thetaMiddle( self )
      res = BiarcMexWrapper( 'thetaMiddle', self.objectHandle );
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
    function [ C0, C1 ] = getCircles( self )
      x0     = self.xBegin();
      y0     = self.yBegin();
      theta0 = self.thetaBegin();
      kappa0 = self.kappa0();
      L      = self.length0();
      C0     = CircleArc( x0, y0, theta0, kappa0, L );
      x0     = self.xMiddle();
      y0     = self.yMiddle();
      theta0 = self.thetaMiddle();
      kappa0 = self.kappa1();
      L      = self.length1();
      C1     = CircleArc( x0, y0, theta0, kappa0, L );
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
    function plot( self, npts, varargin )
      if nargin<2
        npts = 64;
      end
      if nargin>2
        fmt1 = varargin{1};
      else
        fmt1 = {'Color','blue','Linewidth',2};
      end
      if nargin>3
        fmt2 = varargin{2};
      else
        fmt2 = {'Color','red','Linewidth',2};
      end
      [C0,C1] = self.getCircles();
      C0.plot(npts,fmt1);
      hold on;
      C1.plot(npts,fmt2);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotCurvature( self, npts, varargin )
      if nargin < 2
        npts = 1000;
      end
      L = BiarcMexWrapper( 'length', self.objectHandle );
      S = 0:L/npts:L;
      [ ~, ~, ~, kappa ] = BiarcMexWrapper( 'evaluate', self.objectHandle, S );
      plot( S, kappa, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotTheta( self, npts, varargin )
      if nargin < 2
        npts = 1000;
      end
      L = BiarcMexWrapper( 'length', self.objectHandle );
      S = 0:L/npts:L;
      [ ~, ~, theta, ~ ] = BiarcMexWrapper( 'evaluate', self.objectHandle, S );
      plot( S, theta, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotNormal( self, step, len )
      for s=0:step:self.length()
        [ x, y, theta, ~ ] = self.evaluate(s);
        n = [sin(theta),-cos(theta)];
        A = [x,x+len*n(1)];
        B = [y,y+len*n(2)];
        plot(A,B);
      end
    end
  end
end
