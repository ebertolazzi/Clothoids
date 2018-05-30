classdef ClothoidList < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function self = ClothoidList()
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function delete( self )
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function push_back( self, c )
      ClothoidListMexWrapper( 'push_back', self.objectHandle, c.obj_handle() );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
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
      [varargout{1:nargout}] = ClothoidListMexWrapper( 'evaluate', self.objectHandle, s );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function varargout = eval( self, s )
      [varargout{1:nargout}] = ClothoidListMexWrapper( 'eval', self.objectHandle, s );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function varargout = eval_D( self, s )
      [varargout{1:nargout}] = ClothoidListMexWrapper( 'eval_D', self.objectHandle, s );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function varargout = eval_DD( self, s )
      [varargout{1:nargout}] = ClothoidListMexWrapper( 'eval_DD', self.objectHandle, s );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function varargout = eval_DDD( self, s )
      [varargout{1:nargout}] = ClothoidListMexWrapper( 'eval_DDD', self.objectHandle, s );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [X,Y,S,DST] = closestPoint( self, qx, qy )
      [X,Y,S,DST] = ClothoidListMexWrapper( 'closestPoint', self.objectHandle, qx, qy );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [DST,S] = distance( self, varargin )
      % eval the angle of the circle curve at curvilinear abscissa `s`
      [DST,S] = ClothoidListMexWrapper( 'distance', self.objectHandle, varargin{:} );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [X,Y,S,DST] = closestPointBySample( self, qx, qy, ds )
      [X,Y,S,DST] = ClothoidListMexWrapper( 'closestPointBySample', self.objectHandle, qx, qy, ds );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [DST,S] = distanceBySample( self, qx, qy, ds )
      [DST,S] = ClothoidListMexWrapper( 'distanceBySample', self.objectHandle, qx, qy, ds );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function res = getX0( self, k )
      res = ClothoidListMexWrapper( 'getX0', self.objectHandle, k );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function res = getY0( self, k )
      res = ClothoidListMexWrapper( 'getY0', self.objectHandle, k );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function res = getTheta0( self, k )
      res = ClothoidListMexWrapper( 'getTheta0', self.objectHandle, k );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function res = getKappa0( self, k )
      res = ClothoidListMexWrapper( 'getKappa0', self.objectHandle, k );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function res = getDkappa( self, k )
      res = ClothoidListMexWrapper( 'getKappa_D', self.objectHandle, k );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function res = getSmin( self, varargin )
      res = ClothoidListMexWrapper( 'getSmin', self.objectHandle, varargin{:}  );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function res = getSmax( self, varargin )
      res = ClothoidListMexWrapper( 'getSmax', self.objectHandle, varargin{:}  );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function res = length( self, varargin )
      res = ClothoidListMexWrapper( 'length', self.objectHandle, varargin{:} );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
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
      ClothoidListMexWrapper( 'rotate', self.objectHandle, angle, cx, cy );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function translate( self, tx, ty )
      % translate the clothoid curve by `(tx,ty)`
      %
      % Usage:
      %    ref.translate(tx, ty)
      %    
      % On input:
      %    tx, ty: horizontal and vertical translation
      %   
      ClothoidListMexWrapper( 'translate', self.objectHandle, tx, ty );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function changeOrigin( self, newX0, newY0 )
      % move the origin of the clothoid to `(newX0, newY0)` 
      %
      % Usage:
      %    ref.changeOrigin(newX0, newY0)
      %    
      % On input:
      %    newX0, newY0: new coordinates of initial point
      %
      ClothoidListMexWrapper( 'changeOrigin', self.objectHandle, newX0, newY0 );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function scale( self, s )
      % scale clothoid by `sc` factor
      %
      % Usage:
      %    ref.scale(newX0, newY0)
      %    
      % On input:
      %    newX0, newY0: new coordinates of initial point
            
      ClothoidListMexWrapper( 'scale', self.objectHandle, s );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function reverse( self )
      % reverse the orientation of the clothoid curve 
      % Usage:
      %    ref.reverse()
      %
      ClothoidListMexWrapper( 'reverse', self.objectHandle );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function BB = bbox( self, max_angle, max_size, varargin )
      % point at infinity
      % Usage:
      %    ref.reverse()
      %
      BB = ClothoidListMexWrapper( 'bbox', self.objectHandle, max_angle, max_size, varargin{:} );
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function lineH = plot( self, step, varargin )
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
      L     = ClothoidCurveMexWrapper( 'length', self.objectHandle );
      S     = 0:step:L ;
      [X,Y] = ClothoidCurveMexWrapper( 'eval', self.objectHandle, S );
      lineH = plot(X,Y, varargin{:});
    end
  end

  end
end