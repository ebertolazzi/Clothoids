classdef CircleArc < CurveBase
  %% MATLAB class wrapper for the underlying C++ class

  methods
    %> Create a new C++ class instance for the circle arc object
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref = CircleArc() % create empty circle
    %>   ref = CircleArc( x0, y0, theta0, k0, L ) % circle passing from (x0,y0)
    %>                                            % at angle theta0 with curvature k0
    %>                                            % and length L
    %>
    %> ```
    %>
    %> **On input:**
    %>
    %> - `x0`, `y0`: coordinate of initial point
    %> - `theta0`:   orientation of the circle at initial point
    %> - `k0`:       curvature of the circle at initial point
    %> - `L`:        length of curve from initial to final point
    %>
    %> **On output:**
    %>
    %> - ref: reference handle to the object instance
    %>
    function self = CircleArc( varargin )
      self@CurveBase( 'CircleArcMexWrapper', 'CircleArc' );
      self.objectHandle = CircleArcMexWrapper( 'new', varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build the circle from known parameters
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    ref.build( x0, y0, theta0, k0, L )
    %>
    %> ```
    %>
    %> build a circle passing from (x0,y0) at angle theta0 with curvature and length
    %>
    %> **On input:**
    %>
    %> - `x0`, `y0`: coordinate of initial point
    %> - `theta0`:   orientation of the circle at initial point
    %> - `k0`:       curvature of the circle at initial point
    %> - `L`:        length of curve from initial to final point
    %>
    function build( self, x0, y0, theta0, k0, L )
      CircleArcMexWrapper( 'build', self.objectHandle, x0, y0, theta0, k0, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build the circle from known parameters
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    ref.build( x0, y0, theta0, x1, y1 ); % circle passing to [x0,y0] and [x1,y1]
    %>                                         % with angle theta0 at [x0,y0]
    %>    ref.build( p0, theta0, p1 ); % circle passing to p0 and p1 with angle theta0 at p0
    %>
    %> ```
    %>
    %> **On input:**
    %>
    %> - `x0`, `y0`: coordinate of initial point
    %> - `theta0`:   orientation of the circle at initial point
    %> - `k0`:       curvature of the circle at initial point
    %> - `L`:        length of curve from initial to final point
    %> - `p0`:       2D point
    %> - `p1`:       2D point
    %>
    function ok = build_G1( self, varargin )
      ok = CircleArcMexWrapper( 'build_G1', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build the circle arc given 3 points. The point can be alingned
    %> in this case a degenerate straight arc if build.
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    ref.build_3P( x0, y0, x1, y1, x2, y2 )
    %>    ref.build_3P( p0, p1, p2 )
    %> ```
    %>
    function ok = build_3P( self, varargin )
      ok = CircleArcMexWrapper( 'build_3P', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Scale circle by `sc` factor
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    ref.scale( sc );
    %> ```
    %>
    function scale( self, sc )
      CircleArcMexWrapper( 'scale', self.objectHandle, sc );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Change the origin of the circle curve to `s0` and set arc lenght to `L`
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    ref.change_curvilinear_origin( s0, L );
    %> ```
    %>
    function change_curvilinear_origin( self, s0, L )
      CircleArcMexWrapper( 'change_curvilinear_origin', self.objectHandle, s0, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function changeCurvilinearOrigin( self, s0, L )
      self.change_curvilinear_origin( s0, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> return a nurbs representation of the circle arc
    function nurbs = to_nurbs( self )
      nurbs = CircleArcMexWrapper( 'to_nurbs', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the arc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plot( npts );
    %>
    %>   fmt = {'Color','blue','Linewidth',2};
    %>   ref.plot( npts, fmt );
    %>
    %> ```
    %>
    %> - `npts`: number of sampling points for plotting
    %> - `fmt` : format of the arc
    %>
    function plot( self, npts, fmt )
      if nargin<2
        npts = 128;
      end
      if nargin<3
        fmt = {'Color','blue','Linewidth',2};
      end
      L = self.length();
      S = 0:L/npts:L;
      [X,Y] = self.eval(S);
      plot(X,Y,fmt{:});
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Export circle parameters
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    S = ref.export();
    %>
    %> ```
    %>
    %> - `S.x0`, `S.y0`: initial point of the circle arc
    %> - `S.x1`, `S.y1`: final point of the circle arc
    %> - `S.theta0`:     initial angle of the circle arc
    %> - `S.theta1`:     final angle of the circle arc
    %> - `S.kappa`:      curvature of the circle arc
    %> - `S.L`:          length of the clothoid arc
    %>
    function S = export( self )
      S.x0     = self.x_begin();
      S.y0     = self.y_begin();
      S.theta0 = self.theta_begin();
      S.kappa  = self.kappa_begin();
      S.x1     = self.x_end();
      S.y1     = self.y_end();
      S.theta1 = self.theta_end();
      S.L      = self.length();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the polygon of the NURBS for the arc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plotPolygon();
    %>   ref.plotPolygon( 'Color','blue','Linewidth',2 );
    %>
    %> ```
    %>
    %> - `fmt` : format of the arc
    %>
    function plot_polygon( self, varargin )
      arc = self.to_nurbs();
      xx  = arc.coefs(1,:)./arc.coefs(3,:);
      yy  = arc.coefs(2,:)./arc.coefs(3,:);
      if nargin > 1
        plot(xx,yy,varargin{:});
      else
        plot(xx,yy,':ok');
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function plotPolygon( self, varargin )
      self.plot_polygon( varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
