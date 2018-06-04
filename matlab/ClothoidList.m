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
    function obj = obj_handle( self )
      obj = self.objectHandle ;
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
    function push_back_G1( self, varargin )
      ClothoidListMexWrapper( 'push_back_G1', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function C = get( self, k )
      [ x0, y0, theta0, k0, dk, L ] = ClothoidListMexWrapper( 'get', self.objectHandle, k ) ;
      C = ClothoidCurve( x0, y0, theta0, k0, dk, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ s, theta, kappa ] = getSTK( self )
      [ s, theta, kappa ] = ClothoidListMexWrapper( 'getSTK', self.objectHandle ) ;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function N = numSegment( self )
      N = ClothoidListMexWrapper( 'numSegment', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = build3arcG2( self, x0, y0, theta0, kappa0, x1, y1, theta1, kappa1 )
      ok = ClothoidListMexWrapper( 'build3arcG2', self.objectHandle, ...
                                   x0, y0, theta0, kappa0, ...
                                   x1, y1, theta1, kappa1 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = build2arcG2( self, x0, y0, theta0, kappa0, x1, y1, theta1, kappa1 )
      ok = ClothoidListMexWrapper( 'build2arcG2', self.objectHandle, ...
                                    x0, y0, theta0, kappa0, ...
                                    x1, y1, theta1, kappa1 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = buildCLC( self, x0, y0, theta0, kappa0, x1, y1, theta1, kappa1 )
      ok = ClothoidListMexWrapper( 'buildCLC', self.objectHandle, ...
                                   x0, y0, theta0, kappa0, ...
                                   x1, y1, theta1, kappa1 );
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
    function res = xBegin( self, varargin )
      res = ClothoidListMexWrapper( 'xBegin', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = xEnd( self, varargin )
      res = ClothoidListMexWrapper( 'xEnd', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = yBegin( self, varargin )
      res = ClothoidListMexWrapper( 'yBegin', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = yEnd( self, varargin )
      res = ClothoidListMexWrapper( 'yEnd', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = thetaBegin( self, varargin )
      res = ClothoidListMexWrapper( 'thetaBegin', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = thetaEnd( self, varargin )
      res = ClothoidListMexWrapper( 'thetaEnd', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = kappaBegin( self, varargin )
      res = ClothoidListMexWrapper( 'kappaBegin', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = kappaEnd( self, varargin )
      res = ClothoidListMexWrapper( 'kappaEnd', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = kappa_D( self, varargin )
      res = ClothoidListMexWrapper( 'kappa_D', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = sMin( self, varargin )
      res = ClothoidListMexWrapper( 'sMin', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = sMax( self, varargin )
      res = ClothoidListMexWrapper( 'sMax', self.objectHandle, varargin{:} );
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
    function export_table( self, filename )
      ClothoidListMexWrapper( 'export_table', self.objectHandle, filename );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function export_ruby( self, filename )
      ClothoidListMexWrapper( 'export_ruby', self.objectHandle, filename );
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
        npts = varargin{1} ;
      else
        npts = 1000 ;
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
          C.plot( npts, fmt1{:} );
        else
          C.plot( npts, fmt2{:} );
        end
        hold on ;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotCurvature( self, varargin )
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
        npts = varargin{1} ;
      else
        npts = 1000 ;
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
      s0 = 0 ;
      for k=1:self.numSegment()
        C  = self.get(k);
        ss = 0:C.length()/npts:C.length();
        [~,~,~,kappa] = C.evaluate(ss);
        if mod(k,2) == 0
          plot( s0+ss, kappa, fmt1{:} );
        else
          plot( s0+ss, kappa, fmt2{:} );
        end
        s0 = s0+ss(end);
        hold on ;
      end
    end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotAngle( self, varargin )
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
        npts = varargin{1} ;
      else
        npts = 1000 ;
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
      s0 = 0 ;
      for k=1:self.numSegment()
        C  = self.get(k);
        ss = 0:C.length()/npts:C.length();
        [~,~,theta,~] = C.evaluate(ss);
        if mod(k,2) == 0
          plot( s0+ss, theta, fmt1{:} );
        else
          plot( s0+ss, theta, fmt2{:} );
        end
        s0 = s0+ss(end);
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