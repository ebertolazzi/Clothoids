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

    function delete(self)
      %% Destroy the C++ class instance
      LineSegmentMexWrapper('delete', self.objectHandle );
    end
        
    function build( self, varargin )
      %
      % Build the circle from known parameters
      %
      % Usage:
      %    (1) ref.build( x0, y0, theta0, k0, L )
      %    (2) ref.build( x0, y0, theta0, k0, smin, smax )
      %    (3) ref.build( p0, p1, p2 )
      %    (4) ref.build( p0, theta0, p1 )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the circle at initial point
      %    k0:     curvature of the circle at initial point
      %    L:      length of curve from initial to final point
      %    smin:   initial curvilinear coordinate of the curve
      %    smax:   final curvilinear coordinate of the curve
      %    p0:     2D point
      %    p1:     2D point
      %    p2:     2D point
      %
      %    (1) circle passing from (x0,y0) at angle theta0 with curvature and length
      %    (2) circle as in (2) with intial and final curvaturevilinear coordinate respect to (x0,y0)
      %    (3) circle arc passing from 3 points
      %    (4) circle passing to p0 and p1 with angle theta0 at p0
      %
      LineSegmentMexWrapper('build', self.objectHandle, varargin{:} );
    end

    function changeOrigin(self, x0, y0)
      %
      %  change the origin of the circlecurve to s0
      %
      %  Usage:
      %    ref.changeOrigin(s0)
      %    
      %  On input:
      %     s0: curvilinear coordinate of the origin of the new curve
      %
      LineSegmentMexWrapper('changeOrigin', self.objectHandle, x0, y0 );
    end

    function translate(self, tx, ty)
      % move the object by `(tx,ty)`
      LineSegmentMexWrapper('translate', self.objectHandle, tx, ty );
    end

    function trim(self, smin, smax)
      % trim circle curve to the corresponding curvilinear parameters
      LineSegmentMexWrapper('trim', self.objectHandle, smin, smax );
    end
        
    function rotate(self, angle, cx, cy)
      % rotate the circle by angle with center of rotation `(cx,cy)`
      % Usage:
      %    ref.rotate(angle, cx, cy)
      %    
      % On input:
      %    angle: the angle of rotation
      %    cx, cy: coordinates of the centre of rotation
      %
      LineSegmentMexWrapper('rotate', self.objectHandle, angle, cx, cy );
    end

    function moveOrigin(self, newX0, newY0)
      % move the origin of the circle to `(newX0, newY0)` 
      % Usage:
      %    ref.moveOrigin(newX0, newY0)
      %
      % On input:
      %    newX0, newY0: new coordinates of initial point
      %
      LineSegmentMexWrapper('moveOrigin', self.objectHandle, newX0, newY0 );
    end

    function nurbs = to_nurbs(self)
      % return a nurbs representation of the circle arc
      nurbs = LineSegmentMexWrapper('to_nurbs', self.objectHandle);
    end

    %% Eval
    function [X,Y] = eval(self, s)
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
      [X,Y] = LineSegmentMexWrapper('eval', self.objectHandle, s );
    end

    function [DX,DY] = eval_D(self, s)
      % eval the circle derivative at curvilinear abscissa `s`
      [DX,DY] = LineSegmentMexWrapper('eval_D', self.objectHandle, s );
    end

    function [DDX,DDY] = eval_DD(self, s)
      % eval the circle second derivative at curvilinear abscissa `s`
      [DDX,DDY] = LineSegmentMexWrapper('eval_DD', self.objectHandle, s );
    end

    %% Eval
    function [DDDX,DDDY] = eval_DDD(self, s)
      % eval the circle third derivative at curvilinear abscissa `s`
      [DDDX,DDDY] = LineSegmentMexWrapper('eval_DDD', self.objectHandle, s );
    end
 
    function X0 = getX0(self)
      X0 = LineSegmentMexWrapper('getX0', self.objectHandle );
    end
 
    function Y0 = getY0(self)
      Y0 = LineSegmentMexWrapper('getY0', self.objectHandle );
    end
 
    function th0 = getTheta0(self)
      th0 = LineSegmentMexWrapper('getTheta0', self.objectHandle );
    end
 
    function res = length(self)
      res = LineSegmentMexWrapper('length', self.objectHandle );
    end

    %% Utils
    function lineH = plot(self,varargin)
      % plot: method to plot the clothoid curve
      % Usage:
      %    lineH = ref.plot()
      %    lineH = ref.plot(0.2)
      %    lineH = ref.plot(0.2, 'r', 'LineWidth', 3)
      %    
      % On input:
      %    step:     the distance between consecutive samples used to
      %              draw the curve
      %    varargin: optional arguments passed to the MATLAB plot
      %              function
      %
      % On output:
      %    lineH: the handle to the line object 
      %           (i.e. the output of the MATLAB plot function)
      seg = self.to_nurbs() ;
      breaks = fnbrk(seg,'b');
      col = 'b';
      lw  = 1;
      if nargin>1 ; col = varargin{1} ; end
      if nargin>2 ; lw  = varargin{2} ; end
      fnplt(seg,[breaks(1),breaks(end)],col,lw);
      hold on;
      % plot poligono
      xx = seg.coefs(1,:)./seg.coefs(3,:);
      yy = seg.coefs(2,:)./seg.coefs(3,:);
      plot(xx,yy,':ok','Color',col);
    end
  end
end