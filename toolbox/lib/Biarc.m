classdef Biarc < CurveBase
  %> MATLAB class wrapper for the underlying C++ class

  methods
    %> Create a new C++ class instance for the Biarc object
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   self = Biarc();
    %>   self = Biarc( x0, y0, theta0, x1, y1, theta1 );
    %>
    %> ```
    %>
    %> **Optional Arguments:**
    %>
    %> - `x0`, `y0`: coordinate of initial point
    %> - `theta0`    : orientation of the clothoid at initial point
    %> - `x1`, `y1`: coordinate of final point
    %> - `theta1`    : orientation of the clothoid at final point
    %>
    %> **On output:**
    %>
    %> - self: reference handle to the object instance
    %>
    function self = Biarc( varargin )
      self@CurveBase( 'BiarcMexWrapper', 'BiArc' );
      self.objectHandle = BiarcMexWrapper( 'new' );
      if nargin > 0
        ok = BiarcMexWrapper( 'build_G1', self.objectHandle, varargin{:} );
        if ~ok
          error('Biarc constructor failed');
        end
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Build the interpolating G1 biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.build_G1( x0, y0, theta0, x1, y1, theta1 );
    %>
    %> ```
    %>
    %> **On input:**
    %>
    %> - `x0`, `y0`: coordinate of initial point
    %> - `theta0`    : orientation of the clothoid at initial point
    %> - `x1`, `y1`: coordinate of final point
    %> - `theta1`    : orientation of the clothoid at final point
    %>
    function ok = build( self, x0, y0, theta0, x1, y1, theta1 )
      ok = BiarcMexWrapper( 'build_G1', self.objectHandle, ...
                            x0, y0, theta0, x1, y1, theta1 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Build the interpolating biarc by 3 points
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.build_3P( x0, y0, x1, y1, x2, y2 );
    %>   ref.build_3P( [x0, y0], [x1, y1], [x2, y2] );
    %>
    %> ```
    %>
    %> **On input:**
    %>
    %> - `x0`, `y0`: coordinate of initial point
    %> - `x1`, `y1`: coordinate of middle point
    %> - `x2`, `y2`: coordinate of final point
    %>
    %> alternative
    %>
    %>    p0: coordinate of initial point
    %>    p1: coordinate of middle point
    %>    p2: coordinate of final point
    %>
    function ok = build_3P( self, varargin )
      ok = BiarcMexWrapper( 'build_3P', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get junction point of the biarc x-coordinate
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   x = ref.x_middle();
    %>
    %> ```
    %>
    function res = x_middle( self )
      res = BiarcMexWrapper( 'x_middle', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function res = xMiddle( self )
      res = self.x_middle();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get junction point of the biarc y-coordinate
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   x = ref.y_middle();
    %>
    %> ```
    %>
    function res = y_middle( self )
      res = BiarcMexWrapper( 'y_middle', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function res = yMiddle( self )
      res = self.y_middle();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get junction point angle of the biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   theta = ref.theta_middle();
    %>
    %> ```
    %>
    function res = theta_middle( self )
      res = BiarcMexWrapper( 'theta_middle', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function res = thetaMiddle( self )
      res = self.theta_middle();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get curvature of the first arc of the biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   kappa0 = ref.kappa0();
    %>
    %> ```
    %>
    function res = kappa0( self )
      res = BiarcMexWrapper( 'kappa0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get curvature of the second arc of the biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   kappa1 = ref.kappa1();
    %>
    %> ```
    %>
    function res = kappa1( self )
      res = BiarcMexWrapper( 'kappa1', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get length of the first arc of the biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   length0 = ref.length0();
    %>
    %> ```
    %>
    function res = length0( self )
      res = BiarcMexWrapper( 'length0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get length of the second arc of the biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   length1 = ref.length1();
    %>
    %> ```
    %>
    function res = length1( self )
      res = BiarcMexWrapper( 'length1', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get the biarc G1 data
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   [x0,y0,theta0,x1,y1,theta1] = ref.getData();
    %>
    %> ```
    %>
    %> - `x0`, `y0`: initial point of the biarc
    %> - `theta0`    : initial angle of the biarc
    %> - `x1`, `y1`: final point of the biarc
    %> - `theta1`    : final angle of the biarc
    %>
    function [ x0,y0,theta0,x1,y1,theta1] = getData( self )
      x0     = self.x_begin();
      y0     = self.y_begin();
      theta0 = self.theta_begin();
      x1     = self.x_end();
      y1     = self.y_end();
      theta1 = self.theta_end();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get the two circle arc composing a biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   [C0,C1] = ref.get_circles();
    %>
    %> ```
    %>
    function [ C0, C1 ] = get_circles( self )
      x0     = self.x_begin();
      y0     = self.y_begin();
      theta0 = self.theta_begin();
      kappa0 = self.kappa0();
      L      = self.length0();
      C0     = CircleArc( x0, y0, theta0, kappa0, L );
      x0     = self.x_middle();
      y0     = self.y_middle();
      theta0 = self.theta_middle();
      kappa0 = self.kappa1();
      L      = self.length1();
      C1     = CircleArc( x0, y0, theta0, kappa0, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ C0, C1 ] = getCircles( self )
      [ C0, C1 ] = self.get_circles();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Return the nurbs represantation of the two arc composing the biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   [arc0,arc1] = ref.to_nurbs();
    %>
    %> ```
    %>
    %> - `arc0`: the nurbs of the first arc
    %> - `arc1`: the nurbs of the second arc
    %>
    function [ arc0, arc1 ] = to_nurbs( self )
      [ arc0, arc1 ] = BiarcMexWrapper( 'to_nurbs', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get the curvilinear coordinates of the point `(x,y)`
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   [s,t] = ref.find_coord( x, y );
    %>
    %> ```
    %>
    %> - `s`: curvilinear coordinate along the curve
    %> - `t`: curvilinear coordinate along the normal of the curve
    %>
    function [ s, t ] = find_coord( self, x, y )
      [ s, t ] = BiarcMexWrapper( 'findST', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plot( npts );
    %>
    %>   fmt1 = {'Color','blue','Linewidth',2};
    %>   fmt2 = {'Color','red','Linewidth',2};
    %>   ref.plot( npts, fmt1, fmt2 );
    %>
    %> ```
    %>
    %> - `npts`: number of sampling points for plotting
    %> - `fmt1`: format of the first arc
    %> - `fmt2`: format of the second arc
    %>
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
      [C0,C1] = self.get_circles();
      C0.plot(npts,fmt1);
      hold on;
      C1.plot(npts,fmt2);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the curvature of the biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plotCurvature( npts );
    %>
    %>   fmt1 = {'Color','blue','Linewidth',2};
    %>   fmt2 = {'Color','red','Linewidth',2};
    %>   ref.plot_curvature( npts, fmt1, fmt2 );
    %>
    %> ```
    %>
    %> - `npts`: number of sampling points for plotting
    %> - `fmt1`: format of the first arc
    %> - `fmt2`: format of the second arc
    %>
    function plot_curvature( self, npts, varargin )
      if nargin < 2
        npts = 1000;
      end
      L = BiarcMexWrapper( 'length', self.objectHandle );
      S = 0:L/npts:L;
      [ ~, ~, ~, kappa ] = BiarcMexWrapper( 'evaluate', self.objectHandle, S );
      plot( S, kappa, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function plotCurvature( self, npts, varargin )
      self.plot_curvature( npts, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the angle of the biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plot_angle( npts );
    %>
    %>   fmt1 = {'Color','blue','Linewidth',2};
    %>   fmt2 = {'Color','red','Linewidth',2};
    %>   ref.plot_angle( npts, fmt1, fmt2 );
    %>
    %> ```
    %>
    %> - `npts`: number of sampling points for plotting
    %> - `fmt1`: format of the first arc
    %> - `fmt2`: format of the second arc
    %>
    function plot_angle( self, npts, varargin )
      if nargin < 2
        npts = 1000;
      end
      L = BiarcMexWrapper( 'length', self.objectHandle );
      S = 0:L/npts:L;
      [ ~, ~, theta, ~ ] = BiarcMexWrapper( 'evaluate', self.objectHandle, S );
      plot( S, theta, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function plotAngle( self, npts, varargin )
      self.plot_angle( npts, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the normal of the biarc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plot_normal( step, len );
    %>
    %> ```
    %>
    %> - `step`: number of sampling normals
    %> - `len`:  length of the plotted normal
    %>
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot_normal( self, step, len )
      for s=0:step:self.length()
        [ x, y, theta, ~ ] = self.evaluate(s);
        n = [sin(theta),-cos(theta)];
        A = [x,x+len*n(1)];
        B = [y,y+len*n(2)];
        plot(A,B);
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function plotNormal( self, step, len )
      self.plot_normal( step, len );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
