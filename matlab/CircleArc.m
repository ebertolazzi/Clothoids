classdef CircleArc < CurveBase
  %% MATLAB class wrapper for the underlying C++ class

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
      self@CurveBase( 'CircleArcMexWrapper' );
      self.objectHandle = CircleArcMexWrapper( 'new', varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function str = is_type( ~ )
      str = 'CircleArc';
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
      CircleArcMexWrapper( 'build', self.objectHandle, x0, y0, theta0, k0, L );
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
      ok = CircleArcMexWrapper( 'build_3P', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function scale( self, sc )
      % scale circle by `sc` factor
      CircleArcMexWrapper( 'scale', self.objectHandle, sc );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function changeCurvilinearOrigin( self, s0, L )
      % change the origin of the circle curve to `s0`
      CircleArcMexWrapper( 'changeCurvilinearOrigin', self.objectHandle, s0, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [p0,p1,p2] = bbTriangles( self, varargin )
      % return the bounding box triangles of the circle arc
      [p0,p1,p2] = CircleArcMexWrapper( 'bbTriangles', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function nurbs = to_nurbs( self )
      % return a nurbs representation of the circle arc
      nurbs = CircleArcMexWrapper( 'to_nurbs', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
    function plotPolygon( self, varargin )
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
    function plotTriangles( self, varargin )
      [p1,p2,p3] = self.bbTriangles();
      for k=1:size(p1,2)
        x = [ p1(1,k), p2(1,k), p3(1,k), p1(1,k) ];
        y = [ p1(2,k), p2(2,k), p3(2,k), p1(2,k) ];
        fill( x, y, 'red','FaceAlpha', 0.5 );
        %plot( x, y, varargin{:});
      end
    end
  end
end
