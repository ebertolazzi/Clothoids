classdef Biarc < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end
    
  methods
    function this = Biarc( varargin )
      %% Create a new C++ class instance for the clothoid arc object
      % Usage:
      %    ref = Biarc()
      %    ref = Biarc( x0, y0, theta0, k0, dk, L )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the clothoid at initial point
      %    k0:     curvature of the clothoid at initial point
      %    dk:     derivative of curvature respect to arclength
      %    L:      length of curve from initial to final point
      %
      %  On output:
      %    ref: reference handle to the object instance
      this.objectHandle = BiarcMexWrapper( 'new' );
      if nargin > 0
        BiarcMexWrapper( 'build', this.objectHandle, varargin{:} );
      end
    end

    function delete(this)
      %% Destroy the C++ class instance
      BiarcMexWrapper('delete', this.objectHandle );
    end

    function build( this, x0, y0, theta0, x1, y1, theta1 )
      % Build the interpolating G1 clothoid arc
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
      BiarcMexWrapper('build', this.objectHandle, x0, y0, theta0, x1, y1, theta1 );
    end

    function varargout = evaluate( this, s )
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
      [varargout{1:nargout}] = BiarcMexWrapper('evaluate', this.objectHandle, s );
    end

    function [x,y] = eval( this, varargin )
      [x,y] = BiarcMexWrapper('eval', this.objectHandle, varargin{:} );
    end

    %% Eval
    function [x_D,y_D] = eval_D( this, varargin )
      [x_D,y_D] = BiarcMexWrapper('eval_D', this.objectHandle, varargin{:} );
    end

    %% Eval
    function [x_DD,y_DD] = eval_DD( this, varargin )
      [x_DD,y_DD] = BiarcMexWrapper('eval_DD', this.objectHandle, varargin{:} );
    end

    %% Eval
    function [x_DDD,y_DDD] = eval_DDD( this, varargin )
      [x_DDD,y_DDD] = BiarcMexWrapper('eval_DDD', this.objectHandle, varargin{:} );
    end

    %%
    function [X,Y,S,DST] = closestPoint( this, qx, qy )
      [X,Y,S,DST] = BiarcMexWrapper('closestPoint', this.objectHandle, qx, qy );
    end

    function [DST,S] = distance( self, varargin )
      % eval the angle of the circle curve at curvilinear abscissa `s`
      [DST,S] = BiarcMexWrapper('distance', self.objectHandle, varargin{:} );
    end

    %%
    function [X,Y,S,DST] = closestPointBySample( this, qx, qy, ds )
      [X,Y,S,DST] = BiarcMexWrapper('closestPointBySample', this.objectHandle, qx, qy, ds );
    end

    %%
    function [DST,S] = distanceBySample( this, qx, qy, ds )
      [DST,S] = BiarcMexWrapper('distanceBySample', this.objectHandle, qx, qy, ds );
    end

    function res = getX0( this )
      res = BiarcMexWrapper('getX0', this.objectHandle );
    end

    function res = getY0( this )
      res = BiarcMexWrapper('getY0', this.objectHandle );
    end

    function res = getTheta0( this )
      res = BiarcMexWrapper('getTheta0', this.objectHandle );
    end

    function res = getKappa0( this )
      res = BiarcMexWrapper('getKappa0', this.objectHandle );
    end

    function res = getDkappa( this )
      res = BiarcMexWrapper('getKappa_D', this.objectHandle );
    end

    function res = getSmin( this )
      res = BiarcMexWrapper('getSmin', this.objectHandle );
    end

    function res = getSmax( this )
      res = BiarcMexWrapper('getSmax', this.objectHandle );
    end

    function res = length( this )
      res = BiarcMexWrapper('length', this.objectHandle );
    end

    function rotate( this, angle, cx, cy )
      % rotate the clothoid curve by angle respect to the centre `(cx,cy)`
      %
      % Usage:
      %    ref.rotate(angle, cx, cy)
      %    
      % On input:
      %    angle: the angle of rotation
      %    cx, cy: coordinates of the centre of rotation
      %
      BiarcMexWrapper('rotate', this.objectHandle, angle, cx, cy );
    end

    function translate( this, tx, ty )
      % translate the clothoid curve by `(tx,ty)`
      %
      % Usage:
      %    ref.translate(tx, ty)
      %    
      % On input:
      %    tx, ty: horizontal and vertical translation
      %   
      BiarcMexWrapper('translate', this.objectHandle, tx, ty );
    end

    function changeOrigin( this, newX0, newY0 )
      % move the origin of the clothoid to `(newX0, newY0)` 
      %
      % Usage:
      %    ref.changeOrigin(newX0, newY0)
      %    
      % On input:
      %    newX0, newY0: new coordinates of initial point
      %
      BiarcMexWrapper('changeOrigin', this.objectHandle, newX0, newY0 );
    end

    function scale( this, s )
      % scale clothoid by `sc` factor
      %
      % Usage:
      %    ref.scale(newX0, newY0)
      %    
      % On input:
      %    newX0, newY0: new coordinates of initial point
            
      BiarcMexWrapper('scale', this.objectHandle, s );
    end

    function reverse( this )
      % reverse the orientation of the clothoid curve 
      % Usage:
      %    ref.reverse()
      %
      BiarcMexWrapper('reverse', this.objectHandle );
    end

    function [xp,yp,xm,ym] = infinity( this )
      % point at infinity
      % Usage:
      %    ref.reverse()
      %
      [xp,yp,xm,ym] = BiarcMexWrapper('infinity', this.objectHandle );
    end

    %% Utils
    function lineH = plot( this, step, varargin )
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
      if nargin<2
        step = 0.1;
      end
      L     = BiarcMexWrapper('length', this.objectHandle );
      S     = 0:step:L ;
      [X,Y] = BiarcMexWrapper('eval', this.objectHandle, S );
      lineH = plot(X,Y, varargin{:});
    end
  end
end