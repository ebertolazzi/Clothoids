classdef Triangle2D < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    function self = Triangle2D( varargin )
      %% Create a new C++ class instance for the triangle object
      %
      % Usage:
      %    (1) ref = Triangle2D()
      %    (2) ref = Triangle2D( x0, y0, x1, y1, x2, y2 )
      %    (3) ref = Triangle2D( p0, p1, p2 )
      %
      %  On output:
      %    ref: reference handle to the object instance
      %
      self.objectHandle = Triangle2DMexWrapper( 'new', varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
      %% Destroy the C++ class instance
      Triangle2DMexWrapper( 'delete', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function build( self, varargin )
      Triangle2DMexWrapper( 'build', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function translate( self, tx, ty )
      % move the object by `(tx,ty)`
      Triangle2DMexWrapper( 'translate', self.objectHandle, tx, ty );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function scale( self, sc )
      % scale triangle by `sc` factor
      Triangle2DMexWrapper( 'scale', self.objectHandle, sc );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [p1,p2,p3] = points( self )
      [p1,p2,p3] = Triangle2DMexWrapper( 'points', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function rotate( self, angle, cx, cy )
      % rotate the triangle by angle with center of rotation `(cx,cy)`
      % Usage:
      %    ref.rotate(angle, cx, cy)
      %    
      % On input:
      %    angle: the angle of rotation
      %    cx, cy: coordinates of the centre of rotation
      %
      Triangle2DMexWrapper('rotate', self.objectHandle, angle, cx, cy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function icode = isInside( self, x, y )
      icode = Triangle2DMexWrapper( 'isInside', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dst = distanceMin( self, x, y )
      dst = Triangle2DMexWrapper( 'distanceMin', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dst = distanceMax( self, x, y )
      dst = Triangle2DMexWrapper( 'distanceMax', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function info( self )
      Triangle2DMexWrapper( 'info', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, color, varargin )
      [p1,p2,p3] = self.points();
      x = [p1(1),p2(1),p3(1),p1(1)];
      y = [p1(2),p2(2),p3(2),p1(2)];
      plot( x, y, varargin{:} );
      fill( x, y, color );
    end
  end
end