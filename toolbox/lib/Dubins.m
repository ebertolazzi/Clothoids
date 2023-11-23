classdef Dubins < matlab.mixin.Copyable
  %> MATLAB class wrapper for the underlying C++ class

  properties (SetAccess = protected, Hidden = true)
    mexName;
    objectHandle; % Handle to the underlying C++ class instance
    call_delete;
    objectType;
  end

  methods(Access = protected)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Make a deep copy of a curve object
    %>
    %> **Usage**
    %>
    %> \rst
    %> .. code-block:: matlab
    %>
    %>   B = A.copy();
    %>
    %> \endrst
    %>
    %> where `A` is the curve object to be copied.
    %>
    function obj = copyElement( self )
      obj              = copyElement@matlab.mixin.Copyable(self);
      obj.objectHandle = feval( self.mexName, 'make_a_copy', self.objectHandle );
      obj.call_delete  = true;
    end
  end

  methods
    %> Create a new C++ class instance for the Biarc object
    %>
    %> **Usage:**
    %>
    %> \rst
    %> .. code-block:: matlab
    %>
    %>   self = Biarc();
    %>   self = Biarc( x0, y0, theta0, x1, y1, theta1 );
    %>
    %> \endrst
    %>
    %> **Optinal Arguments:**
    %>
    %> - `x0`, `y0`: coordinate of initial point
    %> - `theta0`    : orientation of the clothoid at initial point
    %> - `x1`, `y1`: coordinate of final point
    %> - `theta1`    : orientation of the clothoid at final point
    %>
    %> **On output:**
    %>
    %> - self: reference handle to the object instance
    %>
    function self = Dubins( varargin )
      self.mexName     = 'DubinsMexWrapper';
      self.objectType  = 'Dubins';
      self.call_delete = true;
      self.objectHandle = DubinsMexWrapper( 'new' );
      if nargin > 0
        ok = DubinsMexWrapper( 'build', self.objectHandle, varargin{:} );
        if ~ok
          error('Dubins constructor failed');
        end
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Build the interpolating G1 biarc
    %>
    %> **Usage:**
    %>
    %> \rst
    %> .. code-block:: matlab
    %>
    %>   ref.build_G1( x0, y0, theta0, x1, y1, theta1 );
    %>
    %> \endrst
    %>
    %> **On input:**
    %>
    %> - `x0`, `y0`: coordinate of initial point
    %> - `theta0`    : orientation of the clothoid at initial point
    %> - `x3`, `y3`: coordinate of final point
    %> - `theta3`    : orientation of the clothoid at final point
    %>
    function ok = build( self, x0, y0, theta0, x3, y3, theta3, k_max )
      ok = DubinsMexWrapper( 'build', self.objectHandle, x0, y0, theta0, x3, y3, theta3, k_max );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function S = get_pars( self )
      S = DubinsMexWrapper( 'get_pars', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [C0,C1,C2] = get_circles( self )
      S = DubinsMexWrapper( 'get_pars', self.objectHandle );
      C0 = CircleArc( S.x0, S.y0, S.theta0, S.kappa0, S.L0 );
      C1 = CircleArc( S.x1, S.y1, S.theta1, S.kappa1, S.L1 );
      C2 = CircleArc( S.x2, S.y2, S.theta2, S.kappa2, S.L2 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the biarc
    %>
    %> **Usage:**
    %>
    %> \rst
    %> .. code-block:: matlab
    %>
    %>   ref.plot( npts );
    %>
    %>   fmt1 = {'Color','blue','Linewidth',2};
    %>   fmt2 = {'Color','red','Linewidth',2};
    %>   ref.plot( npts, fmt1, fmt2 );
    %>
    %> \endrst
    %>
    %> - `npts`: number of sampling points for plotting
    %> - `fmt1`: format of the first arc
    %> - `fmt2`: format of the second arc
    %>
    function plot( self, npts, varargin )
      if nargin<2
        npts = 64;
      end
      if nargin>2
        fmt1 = varargin{1};
      else
        fmt1 = {'Color','blue','Linewidth',2};
      end
      if nargin>3
        fmt2 = varargin{2};
      else
        fmt2 = {'Color','red','Linewidth',2};
      end
      [C0,C1,C2] = self.get_circles();
      C0.plot(npts,fmt1);
      hold on;
      C1.plot(npts,fmt2);
      C2.plot(npts,fmt1);
    end
  end
end
