 classdef Dubins3p < CurveBase
  %> MATLAB class wrapper for the underlying C++ class

  %properties (SetAccess = protected, Hidden = true)
  %  mexName;
  %  objectHandle; % Handle to the underlying C++ class instance
  %  call_delete;
  %  objectType;
  %end

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
    %function obj = copyElement( self )
    %  obj              = copyElement@matlab.mixin.Copyable(self);
    %  obj.objectHandle = feval( self.mexName, 'make_a_copy', self.objectHandle );
    %  obj.call_delete  = true;
    %end
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
    function self = Dubins3p( varargin )
      self@CurveBase( 'Dubins3pMexWrapper', 'Dubins' );
      self.objectHandle = Dubins3pMexWrapper( 'new' );
      if nargin > 0
        ok = Dubins3pMexWrapper( 'build', self.objectHandle, varargin{:} );
        if ~ok
          error('Dubins3p constructor failed');
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
    function ok = build( self, x0, y0, theta0, xm, ym, x1, y1, theta1, k_max, method )
      ok = Dubins3pMexWrapper( 'build', self.objectHandle, x0, y0, theta0, xm, ym, x1, y1, theta1, k_max, method );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function S = get_pars( self )
      S = Dubins3pMexWrapper( 'get_pars', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ C1, C2, C3, C4, C5, C6 ] = get_circles( self )
      S  = Dubins3pMexWrapper( 'get_pars', self.objectHandle );
      C1 = CircleArc( S.x0, S.y0, S.theta0, S.kappa1, S.L1 );
      C2 = CircleArc( S.x1, S.y1, S.theta1, S.kappa2, S.L2 );
      C3 = CircleArc( S.x2, S.y2, S.theta2, S.kappa3, S.L3 );
      C4 = CircleArc( S.x3, S.y3, S.theta3, S.kappa4, S.L4 );
      C5 = CircleArc( S.x4, S.y4, S.theta4, S.kappa5, S.L5 );
      C6 = CircleArc( S.x5, S.y5, S.theta5, S.kappa6, S.L6 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [L,R] = curve_type( self )
      [L,R] = Dubins3pMexWrapper( 'curve_type', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [l,s] = curve_type_string( self )
      % l = long string
      % s = short string
      [l,s] = Dubins3pMexWrapper( 'curve_type_string', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = num_evaluation( self )
      res = Dubins3pMexWrapper( 'num_evaluation', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function set_tolerance( self, tol )
      Dubins3pMexWrapper( 'set_tolerance', self.objectHandle, tol );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function set_max_evaluation( self, max_eval )
      Dubins3pMexWrapper( 'set_max_evaluation', self.objectHandle, max_eval );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function set_sample_angle( self, ang )
      Dubins3pMexWrapper( 'set_sample_angle', self.objectHandle, ang );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function sample_points( self, npts )
      Dubins3pMexWrapper( 'sample_points', self.objectHandle, npts );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the Dubins
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
      [ C0, C1, C2, C3, C4, C5 ] = self.get_circles();
      C0.plot(npts,fmt1);
      hold on;
      C1.plot(npts,fmt2);
      C2.plot(npts,fmt1);
      C3.plot(npts,fmt2);
      C4.plot(npts,fmt1);
      C5.plot(npts,fmt2);
    end
  end
end
