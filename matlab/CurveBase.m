classdef CurveBase < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = protected, Hidden = true)
    mexName;
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    function self = CurveBase( mexName )
      self.mexName = mexName;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
      %% Destroy the C++ class instance
      if self.objectHandle ~= 0
        feval( self.mexName, 'delete', self.objectHandle );
      end
      self.objectHandle = 0; % avoid double destruction of object
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function obj = obj_handle( self )
      obj = self.objectHandle;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function copy( self, C )
      feval( self.mexName, 'copy', self.objectHandle, C.obj_handle() );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ xmin, ymin, xmax, ymax ] = bbox( self, varargin )
      % return the bounding box triangle of the circle arc
      [ xmin, ymin, xmax, ymax ] = ...
        feval( self.mexName, 'bbox', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function translate( self, tx, ty )
      % translate curve by vector (tx,ty)
      feval( self.mexName, 'translate', self.objectHandle, tx, ty );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function trim( self, smin, smax )
      % trim circle curve to the corresponding curvilinear parameters
      feval( self.mexName, 'trim', self.objectHandle, smin, smax );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function rotate( self, angle, cx, cy )
      feval( self.mexName, 'rotate', self.objectHandle, angle, cx, cy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function reverse( self )
      feval( self.mexName, 'reverse', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function scale( self, sc )
      feval( self.mexName, 'scale', self.objectHandle, sc );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function changeOrigin( self, newX0, newY0 )
      feval( self.mexName, 'changeOrigin', self.objectHandle, newX0, newY0 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x, y, theta, kappa ] = evaluate( self, s, varargin )
      [ x, y, theta, kappa ] = feval( self.mexName, 'evaluate', self.objectHandle, s, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval( self, varargin )
      [ varargout{1:nargout} ] = ...
        feval( self.mexName, 'eval', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_D( self, varargin )
      % eval the circle derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        feval( self.mexName, 'eval_D', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DD( self, varargin )
      % eval the circle second derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        feval( self.mexName, 'eval_DD', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DDD( self, varargin )
      % eval the circle third derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        feval( self.mexName, 'eval_DDD', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = theta(self, s)
      th = feval( self.mexName, 'theta', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = theta_D(self, s)
      th = feval( self.mexName, 'theta_D', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = theta_DD(self, s)
      th = feval( self.mexName, 'theta_DD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = theta_DDD(self, s)
      th = feval( self.mexName, 'theta_DDD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = kappa(self, s)
      th = feval( self.mexName, 'kappa', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = kappa_D(self, s)
      th = feval( self.mexName, 'kappa_D', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th = kappa_DD(self, s)
      th = feval( self.mexName, 'kappa_DD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x, y ] = xyBegin( self )
      [ x, y ] = feval( self.mexName, 'xyBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function X0 = xBegin( self )
      X0 = feval( self.mexName, 'xBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Y0 = yBegin( self )
      Y0 = feval( self.mexName, 'yBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th0 = thetaBegin( self )
      th0 = feval( self.mexName, 'thetaBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th0 = kappaBegin( self )
      th0 = feval( self.mexName, 'kappaBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x, y ] = xyEnd( self )
      [ x, y ] = feval( self.mexName, 'xyEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function X1 = xEnd( self )
      X1 = feval( self.mexName, 'xEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Y1 = yEnd( self )
      Y1 = feval( self.mexName, 'yEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th1 = thetaEnd( self )
      th1 = feval( self.mexName, 'thetaEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th1 = kappaEnd( self )
      th1 = feval( self.mexName, 'kappaEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = length( self, varargin )
      res = feval( self.mexName, 'length', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ p1, p2 ] = points( self )
      [ p1, p2 ] = feval( self.mexName, 'points', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x, y, s, t, iflag, dst ] = closestPoint( self, qx, qy, varargin )
      [ x, y, s, t, iflag, dst ] = ...
        feval( self.mexName, 'closestPoint', self.objectHandle, qx, qy, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dst = distance( self, qx, qy, varargin )
      dst = feval( self.mexName, 'distance', self.objectHandle, qx, qy, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = collision( self, OBJ, varargin )
      ok = feval( self.mexName, 'collision', ...
                                self.objectHandle, ...
                                OBJ.obj_handle(), OBJ.is_type(), ...
                                varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [s1,s2] = intersect( self, OBJ, varargin )
      [s1,s2] = feval( self.mexName, 'intersect', ...
                                     self.objectHandle, ...
                                     OBJ.obj_handle(), OBJ.is_type(), ...
                                     varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function info( self )
      feval( self.mexName, 'info', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [s,t] = find_coord( self, x, y )
      [s,t] = feval( self.mexName, 'findST', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function yesAABBtree( self )
      feval( self.mexName, 'yesAABBtree', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function noAABBtree( self )
      feval( self.mexName, 'noAABBtree', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotTBox( self, P1, P2, P3, varargin )
      for k=1:size(P1,2)
        pp1 = P1(:,k);
        pp2 = P2(:,k);
        pp3 = P3(:,k);
        plot( [pp1(1),pp2(1),pp3(1),pp1(1)], ...
              [pp1(2),pp2(2),pp3(2),pp1(2)], varargin{:} );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotBBox( self, varargin )
      if nargin > 1
        offs = varargin{1};
        [xmin,ymin,xmax,ymax] = self.bbox( offs );
      else
        [xmin,ymin,xmax,ymax] = self.bbox();
      end
      x = [ xmin, xmax, xmax, xmin, xmin ];
      y = [ ymin, ymin, ymax, ymax, ymin ];
      if nargin > 2
        plot( x, y, varargin{2:end});
      else
        plot( x, y );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
