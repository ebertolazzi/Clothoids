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
      self.objectHandle = CircleMexWrapper('new', varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      %% Destroy the C++ class instance
      CircleMexWrapper('delete', self.objectHandle );
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
      CircleMexWrapper('build', self.objectHandle, x0, y0, theta0, k0, L );
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
      ok = CircleMexWrapper('build_G1', self.objectHandle, varargin{:} );
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
      ok = CircleMexWrapper('build_3P', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function translate(self, tx, ty)
      % move the object by `(tx,ty)`
      CircleMexWrapper('translate', self.objectHandle, tx, ty );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function changeCurvilinearOrigin(self, s0, L )
      % change the origin of the circle curve to `s0`
      CircleMexWrapper('changeCurvilinearOrigin', self.objectHandle, s0, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function trim(self, smin, smax)
      % trim circle curve to the corresponding curvilinear parameters
      CircleMexWrapper('trim', self.objectHandle, smin, smax );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function scale(self, sc)
      % scale circle by `sc` factor
      CircleMexWrapper('scale', self.objectHandle, sc );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function rotate(self, angle, cx, cy)
      % rotate the circle by angle with center of rotation `(cx,cy)`
      % Usage:
      %    ref.rotate(angle, cx, cy)
      %    
      % On input:
      %    angle: the angle of rotation
      %    cx, cy: coordinates of the centre of rotation
      %
      CircleMexWrapper('rotate', self.objectHandle, angle, cx, cy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function changeOrigin(self, newX0, newY0)
      % move the origin of the circle to `(newX0, newY0)` 
      % Usage:
      %    ref.changeOrigin(newX0, newY0)
      %
      % On input:
      %    newX0, newY0: new coordinates of initial point
      %
      CircleMexWrapper('changeOrigin', self.objectHandle, newX0, newY0 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [p0,p1,p2,ok] = bbTriangle(self)
      % return the bounding box triangle of the circle arc
      [p0,p1,p2,ok] = CircleMexWrapper('bbTriangle', self.objectHandle);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function nurbs = to_nurbs(self)
      % return a nurbs representation of the circle arc
      nurbs = CircleMexWrapper('to_nurbs', self.objectHandle);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
      [X,Y] = CircleMexWrapper('eval', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [DX,DY] = eval_D(self, s)
      % eval the circle derivative at curvilinear abscissa `s`
      [DX,DY] = CircleMexWrapper('eval_D', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [DDX,DDY] = eval_DD(self, s)
      % eval the circle second derivative at curvilinear abscissa `s`
      [DDX,DDY] = CircleMexWrapper('eval_DD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [DDDX,DDDY] = eval_DDD(self, s)
      % eval the circle third derivative at curvilinear abscissa `s`
      [DDDX,DDDY] = CircleMexWrapper('eval_DDD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = theta(self, s)
      % eval the angle of the circle curve at curvilinear abscissa `s`
      th = CircleMexWrapper('theta', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function X0 = xBegin(self)
      X0 = CircleMexWrapper('xBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function X1 = xEnd(self)
      X1 = CircleMexWrapper('xEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Y0 = yBegin(self)
      Y0 = CircleMexWrapper('yBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Y1 = yEnd(self)
      Y1 = CircleMexWrapper('yEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th0 = thetaBegin(self)
      th0 = CircleMexWrapper('thetaBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th1 = thetaEnd(self)
      th1 = CircleMexWrapper('thetaEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function kappa = getKappa(self)
      kappa = CircleMexWrapper('getKappa', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = length(self)
      res = CircleMexWrapper('length', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [DST,S] = distance( self, varargin )
      % eval the angle of the circle curve at curvilinear abscissa `s`
      [DST,S] = CircleMexWrapper('distance', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function lineH = plot(self,varargin)
      % plot: method to plot the circle curve
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
      arc = self.to_nurbs() ;
      breaks = fnbrk(arc,'b');
      col = 'b';
      lw  = 1;
      if nargin>1 ; col = varargin{1} ; end
      if nargin>2 ; lw  = varargin{2} ; end
      fnplt(arc,[breaks(1),breaks(end)],col,lw);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function lineH = plotPolygon(self,varargin)
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