classdef Triangle2D < matlab.mixin.Copyable
  %> MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
    call_delete;
  end

  methods
    %>
    %> Create a new C++ class instance for the triangle object
    %>
    %> **Usage:**
    %>
    %> \rst
    %> .. code-block:: matlab
    %>
    %>    ref = Triangle2D();
    %>    ref = Triangle2D( x0, y0, x1, y1, x2, y2 );
    %>    ref = Triangle2D( p0, p1, p2 );
    %>
    %> \endrst
    %>
    %> **On output:**
    %>
    %> - `ref`: reference handle to the object instance
    %>
    function self = Triangle2D( varargin )
      self.objectHandle = Triangle2DMexWrapper( 'new', varargin{:} );
      self.call_delete  = true;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Destroy the C++ class instance
    function delete( self )
      %% Destroy the C++ class instance
      if self.objectHandle ~= 0
        if self.call_delete
          Triangle2DMexWrapper( 'delete', self.objectHandle );
          self.objectHandle = 0; % avoid double destruction of object
        end
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function obj = obj_handle( self )
      obj = self.objectHandle;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function str = is_type( ~ )
      str = 'Triangle2D';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function build( self, varargin )
      Triangle2DMexWrapper( 'build', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> move the triangle by `(tx,ty)`
    function translate( self, tx, ty )
      Triangle2DMexWrapper( 'translate', self.objectHandle, tx, ty );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> scale triangle by `sc` factor
    function scale( self, sc )
      Triangle2DMexWrapper( 'scale', self.objectHandle, sc );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [p1,p2,p3] = points( self )
      [p1,p2,p3] = Triangle2DMexWrapper( 'points', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Rotate the triangle by angle with center of rotation `(cx,cy)`
    %>
    %> **Usage:**
    %>
    %> \rst
    %> .. code-block:: matlab
    %>
    %>    ref.rotate(angle, cx, cy)
    %>
    %> \endrst
    %>
    %> **On input:**
    %>
    %> - `angle`: the angle of rotation
    %> - `cx`, `cy`: coordinates of the centre of rotation
    %>
    function rotate( self, angle, cx, cy )
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
    function yesno = overlap( self, obj )
      yesno = Triangle2DMexWrapper( 'overlap', self.objectHandle, obj.obj_handle() );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, color, varargin )
      [p1,p2,p3] = self.points();
      x = [p1(1),p2(1),p3(1),p1(1)];
      y = [p1(2),p2(2),p3(2),p1(2)];
      fill( x, y, color );
      plot( x, y, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot2( self, varargin )
      [p1,p2,p3] = self.points();
      x = [p1(1),p2(1),p3(1),p1(1)];
      y = [p1(2),p2(2),p3(2),p1(2)];
      fill( x, y, varargin{:} );
    end
  end
end
