classdef PolyLine < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    function self = PolyLine( )
      self.objectHandle = PolyLineMexWrapper( 'new' );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      %% Destroy the C++ class instance
      PolyLineMexWrapper( 'delete', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function obj = obj_handle( self )
      obj = self.objectHandle ;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function str = is_type( ~ )
      str = 'PolyLine' ;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function copy( C )
      PolyLineMexWrapper( 'copy', self.objectHandle, C.obj_handle() );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function build( self, x, y )
      PolyLineMexWrapper( 'build', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function approx( self, obj, tol )
      PolyLineMexWrapper( 'approx', self.objectHandle, ...
                                    obj.objectHandle,  ...
                                    tol,               ...
                                    obj.is_type() );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [s1,s2] = intersect( self, obj )
      [s1,s2] = PolyLineMexWrapper( 'intersect', ...
                                    self.objectHandle, ...
                                    obj.objectHandle,  ...
                                    obj.is_type() );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function translate( self, tx, ty )
      % move the object by `(tx,ty)`
      PolyLineMexWrapper( 'translate', self.objectHandle, tx, ty );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function rotate( self, angle, cx, cy )
      % rotate the circle by angle with center of rotation `(cx,cy)`
      % Usage:
      %    ref.rotate(angle, cx, cy)
      %
      % On input:
      %    angle: the angle of rotation
      %    cx, cy: coordinates of the centre of rotation
      %
      PolyLineMexWrapper( 'rotate', self.objectHandle, angle, cx, cy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function reverse( self )
      PolyLineMexWrapper( 'reverse', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval( self, s )
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
      [ varargout{1:nargout} ] = ...
        PolyLineMexWrapper( 'eval', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_D( self, s )
      % eval the circle derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        PolyLineMexWrapper( 'eval_D', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DD( self, s )
      % eval the circle second derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        PolyLineMexWrapper( 'eval_DD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DDD( self, s )
      % eval the circle third derivative at curvilinear abscissa `s`
      [ varargout{1:nargout} ] = ...
        PolyLineMexWrapper( 'eval_DDD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function X0 = xBegin( self )
      X0 = PolyLineMexWrapper( 'xBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Y0 = yBegin( self )
      Y0 = PolyLineMexWrapper( 'yBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function th0 = theta( self )
      th0 = PolyLineMexWrapper( 'theta', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = length( self )
      res = PolyLineMexWrapper( 'length', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x, y ] = polygon( self )
      [ x, y ] = PolyLineMexWrapper( 'polygon', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ d, s ] = distance( self, x, y )
      [ d, s ] = PolyLineMexWrapper( 'distance', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, varargin )
      if nargin > 1
        fmt1 = varargin{1} ;
      else
        fmt1 = {'Color','red','LineWidth',3} ;
      end
      if nargin > 2
        fmt2 = varargin{2} ;
      else
        fmt2 = {'Color','blue','LineWidth',3} ;
      end
      [ x, y ] = self.polygon();
      for k=2:length(x)
        if mod(k,2) == 0
          plot( x(k-1:k), y(k-1:k), fmt1{:} );
        else
          plot( x(k-1:k), y(k-1:k), fmt2{:} );
        end
        hold on;
      end
    end
  end
end
