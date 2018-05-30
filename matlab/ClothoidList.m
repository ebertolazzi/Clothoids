classdef ClothoidList < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = ClothoidList()
      self.objectHandle = ClothoidListMexWrapper( 'new' );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function reserve( self, N )
      ClothoidListMexWrapper( 'reserve', self.objectHandle, N );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function push_back( self, varargin )
      if nargin == 2
        ClothoidListMexWrapper( 'push_back', self.objectHandle, varargin{1}.obj_handle() );
      else
        ClothoidListMexWrapper( 'push_back', self.objectHandle, varargin{:} );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function C = get( self, k )
      [ x0, y0, theta0, k0, dk, L ] = ClothoidListMexWrapper( 'get', self.objectHandle, k ) ;
      C = ClothoidCurve( x0, y0, theta0, k0, dk, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function N = numSegment( self )
      N = ClothoidListMexWrapper( 'numSegment', self.objectHandle );
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
      [varargout{1:nargout}] = ClothoidListMexWrapper( 'evaluate', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval( self, s )
      [varargout{1:nargout}] = ClothoidListMexWrapper( 'eval', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_D( self, s )
      [varargout{1:nargout}] = ClothoidListMexWrapper( 'eval_D', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DD( self, s )
      [varargout{1:nargout}] = ClothoidListMexWrapper( 'eval_DD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DDD( self, s )
      [varargout{1:nargout}] = ClothoidListMexWrapper( 'eval_DDD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [X,Y,S,DST] = closestPoint( self, qx, qy )
      [X,Y,S,DST] = ClothoidListMexWrapper( 'closestPoint', self.objectHandle, qx, qy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [DST,S] = distance( self, varargin )
      % eval the angle of the circle curve at curvilinear abscissa `s`
      [DST,S] = ClothoidListMexWrapper( 'distance', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [X,Y,S,DST] = closestPointBySample( self, qx, qy, ds )
      [X,Y,S,DST] = ClothoidListMexWrapper( 'closestPointBySample', self.objectHandle, qx, qy, ds );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [DST,S] = distanceBySample( self, qx, qy, ds )
      [DST,S] = ClothoidListMexWrapper( 'distanceBySample', self.objectHandle, qx, qy, ds );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = getX0( self, k )
      res = ClothoidListMexWrapper( 'getX0', self.objectHandle, k );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = getY0( self, k )
      res = ClothoidListMexWrapper( 'getY0', self.objectHandle, k );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = getTheta0( self, k )
      res = ClothoidListMexWrapper( 'getTheta0', self.objectHandle, k );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = getKappa0( self, k )
      res = ClothoidListMexWrapper( 'getKappa0', self.objectHandle, k );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = getDkappa( self, k )
      res = ClothoidListMexWrapper( 'getKappa_D', self.objectHandle, k );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = getSmin( self, varargin )
      res = ClothoidListMexWrapper( 'getSmin', self.objectHandle, varargin{:}  );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = getSmax( self, varargin )
      res = ClothoidListMexWrapper( 'getSmax', self.objectHandle, varargin{:}  );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = length( self, varargin )
      res = ClothoidListMexWrapper( 'length', self.objectHandle, varargin{:} );
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
      ClothoidListMexWrapper( 'rotate', self.objectHandle, angle, cx, cy );
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
      ClothoidListMexWrapper( 'translate', self.objectHandle, tx, ty );
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
      ClothoidListMexWrapper( 'changeOrigin', self.objectHandle, newX0, newY0 );
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
            
      ClothoidListMexWrapper( 'scale', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function reverse( self )
      % reverse the orientation of the clothoid curve 
      % Usage:
      %    ref.reverse()
      %
      ClothoidListMexWrapper( 'reverse', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function BB = bbox( self, max_angle, max_size, varargin )
      % point at infinity
      % Usage:
      %    ref.reverse()
      %
      BB = ClothoidListMexWrapper( 'bbox', self.objectHandle, max_angle, max_size, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, varargin )
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
      if nargin > 1
        step = varargin{1} ;
      else
        step = 1000 ;
      end
      if nargin > 2
        fmt1 = varargin{2} ;
      else
        fmt1 = {'Color','red','LineWidth',3} ;
      end
      if nargin > 3
        fmt2 = varargin{3} ;
      else
        fmt2 = {'Color','blue','LineWidth',3} ;
      end
      for k=1:self.numSegment()
        C = self.get(k);
        if mod(k,2) == 0
          C.plot( step, fmt1{:} );
        else
          C.plot( step, fmt2{:} );
        end
        hold on ;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotNormal( self, step, len )
      for k=1:self.numSegment()
        C = self.get(k);
        C.plotNormal( step, len );
      end
    end
  end
  
end