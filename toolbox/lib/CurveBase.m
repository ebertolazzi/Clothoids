classdef CurveBase < matlab.mixin.Copyable
  %% MATLAB class wrapper for the underlying C++ class
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
    %> ```{matlab}
    %>
    %>   B = A.copy();
    %>
    %> ```
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
    function self = CurveBase( mexName, objectType )
      self.mexName     = mexName;
      self.objectType  = objectType;
      self.call_delete = true;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
      %% Destroy the C++ class instance
      if self.objectHandle ~= 0
        if self.call_delete
          feval( self.mexName, 'delete', self.objectHandle );
          self.objectHandle = 0; % avoid double destruction of object
        end
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Return the `pointer` of the interbal stored c++ object
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   obj = ref.obj_handle();
    %>
    %> ```
    %>
    function obj = obj_handle( self )
      obj = self.objectHandle;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function str = is_type( self )
      str = self.objectType;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function load( self, OBJ )
      feval( self.mexName, ...
        'build2', self.objectHandle, OBJ.obj_handle(), OBJ.is_type() ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Return the bounding box of the curve object
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   [ xmin, ymin, xmax, ymax ] = ref.bbox( C );
    %>   [ xmin, ymin, xmax, ymax ] = ref.bbox( C, offs );
    %>   [ xmin, ymin, xmax, ymax ] = ref.bbox( C, offs, 'ISO' );
    %>   [ xmin, ymin, xmax, ymax ] = ref.bbox( C, offs, 'SAE' );
    %>
    %> ```
    %>
    %> - xmin: x minimum coordinate of the bounding box
    %> - ymin: y minimum coordinate of the bounding box
    %> - xmax: x maximum coordinate of the bounding box
    %> - ymax: y maximum coordinate of the bounding box
    %>
    %> **Optional Arguments**
    %>
    %> - offs: offset of the curve used in the bbox computation
    %> - 'ISO'/'SAE': use ISO or SAE orientation of the normal for offset computation
    %>
    function [ xmin, ymin, xmax, ymax ] = bbox( self, varargin )
      % return the bounding box triangle of the circle arc
      [ xmin, ymin, xmax, ymax ] = ...
        feval( self.mexName, 'bbox', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Translate the curve by `(tx,ty)`
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   ref.translate( tx, ty );
    %>
    %> ```
    %>
    function translate( self, tx, ty )
      % translate curve by vector (tx,ty)
      feval( self.mexName, 'translate', self.objectHandle, tx, ty );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Cut the curve at the curvilinear parameter `smin` up to `smax`
    %>
    %> **Usage**
    %>
    %>    ref.trim( smin, smax );
    %>
    function trim( self, smin, smax )
      % trim circle curve to the corresponding curvilinear parameters
      feval( self.mexName, 'trim', self.objectHandle, smin, smax );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Rotate the curve by angle `angle` around point `(cx, cy)`
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   ref.rotate( angle, cx, cy );
    %>
    %> ```
    %>
    function rotate( self, angle, cx, cy )
      % rotate curve around `(cx,cy)` by angle `angle`
      feval( self.mexName, 'rotate', self.objectHandle, angle, cx, cy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Reverse the direction of travel of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   ref.reverse();
    %>
    %> ```
    %>
    function reverse( self )
      % reverse curve mileage
      feval( self.mexName, 'reverse', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Scale the curve by factor `sc`
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   ref.scale( sc );
    %>
    %> ```
    %>
    function scale( self, sc )
      % scale curve by `sc`
      feval( self.mexName, 'scale', self.objectHandle, sc );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Translate the curve in such a way the origin is at `(newX0,newY0)`.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   ref.change_origin( newX0, newY0 );
    %>
    %> ```
    %>
    function change_origin( self, newX0, newY0 )
      % change the origgin or the curve to `(newX0,newY0)`
      feval( self.mexName, 'change_origin', self.objectHandle, newX0, newY0 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function changeOrigin( self, newX0, newY0 )
      self.change_origin(newX0, newY0);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   [ x, y, theta, kappa ] = ref.evaluate( s );
    %>   [ x, y, theta, kappa ] = ref.evaluate( s, offs );
    %>   [ x, y, theta, kappa ] = ref.evaluate( s, offs, 'ISO' );
    %>   [ x, y, theta, kappa ] = ref.evaluate( s, offs, 'SAE' );
    %>
    %> ```
    %>
    %> **Optional Arguments**
    %>
    %> - offs: offset of the curve used in computation
    %> - 'ISO'/'SAE': use ISO or SAE orientation of the normal for offset computation
    %>
    function [ x, y, theta, kappa ] = evaluate( self, s, varargin )
      [ x, y, theta, kappa ] = feval( self.mexName, 'evaluate', ...
         self.objectHandle, s, varargin{:} ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   XY = ref.eval( s );
    %>   XY = ref.eval( s, offs );
    %>   XY = ref.eval( s, offs, 'ISO' );
    %>   XY = ref.eval( s, offs, 'SAE' );
    %>
    %>   [X,Y] = ref.eval( s );
    %>   [X,Y] = ref.eval( s, offs );
    %>   [X,Y] = ref.eval( s, offs, 'ISO' );
    %>   [X,Y] = ref.eval( s, offs, 'SAE' );
    %>
    %> ```
    %>
    %> **Optional Arguments**
    %>
    %> - offs: offset of the curve used compiutation
    %> - 'ISO'/'SAE': use ISO or SAE orientation of the normal for offset computation
    %>
    %> **Output**
    %>
    %> - XY: matrix `2 x n` of the evaluated points
    %> - X: vector of the x-coordinates of the evaluated points
    %> - Y: vector of the y-coordinates of the evaluated points
    %>
    function varargout = eval( self, varargin )
      [ varargout{1:nargout} ] = ...
        feval( self.mexName, 'eval', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the first derivatives of the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   XY = ref.eval_D( s );
    %>   XY = ref.eval_D( s, offs );
    %>   XY = ref.eval_D( s, offs, 'ISO' );
    %>   XY = ref.eval_D( s, offs, 'SAE' );
    %>
    %>   [X,Y] = ref.eval_D( s );
    %>   [X,Y] = ref.eval_D( s, offs );
    %>   [X,Y] = ref.eval_D( s, offs, 'ISO' );
    %>   [X,Y] = ref.eval_D( s, offs, 'SAE' );
    %>
    %> ```
    %>
    %> **Optional Arguments**
    %>
    %> - `offs`: offset of the curve used in computation
    %> - 'ISO'/'SAE': use ISO or SAE orientation of the normal for offset computation
    %>
    %> **Output**
    %>
    %> - `XY`: matrix `2 x n` of the evaluated points
    %> - `X`: vector of the x-coordinates of the evaluated point derivatives
    %> - `Y`: vector of the y-coordinates of the evaluated point derivatives
    %>
    function varargout = eval_D( self, varargin )
      % evaluate derivative at `s`
      % XY = eval_D(s,[offs,'ISO'/'SAE')
      % [X,Y] = eval_D(s,[offs,'ISO'/'SAE')
      [ varargout{1:nargout} ] = ...
        feval( self.mexName, 'eval_D', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the second derivatives of the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   XY = ref.eval_DD( s );
    %>   XY = ref.eval_DD( s, offs );
    %>   XY = ref.eval_DD( s, offs, 'ISO' );
    %>   XY = ref.eval_DD( s, offs, 'SAE' );
    %>
    %>   [X,Y] = ref.eval_DD( s );
    %>   [X,Y] = ref.eval_DD( s, offs );
    %>   [X,Y] = ref.eval_DD( s, offs, 'ISO' );
    %>   [X,Y] = ref.eval_DD( s, offs, 'SAE' );
    %>
    %> ```
    %>
    %> **Optional Arguments**
    %>
    %> - `offs`: offset of the curve used in computation
    %> - 'ISO'/'SAE': use ISO or SAE orientation of the normal for offset computation
    %>
    %> **Output**
    %>
    %> - `XY`: matrix `2 x n` of the evaluated points
    %> - `X`: vector of the x-coordinates of the evaluated point derivatives
    %> - `Y`: vector of the y-coordinates of the evaluated point derivatives
    %>
    function varargout = eval_DD( self, varargin )
      [ varargout{1:nargout} ] = ...
        feval( self.mexName, 'eval_DD', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the third derivatives of the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   XY = ref.eval_DDD( s );
    %>   XY = ref.eval_DDD( s, offs );
    %>   XY = ref.eval_DDD( s, offs, 'ISO' );
    %>   XY = ref.eval_DDD( s, offs, 'SAE' );
    %>
    %>   [X,Y] = ref.eval_DDD( s );
    %>   [X,Y] = ref.eval_DDD( s, offs );
    %>   [X,Y] = ref.eval_DDD( s, offs, 'ISO' );
    %>   [X,Y] = ref.eval_DDD( s, offs, 'SAE' );
    %>
    %> ```
    %>
    %> **Optional Arguments**
    %>
    %> - `offs`: offset of the curve used in computation
    %> - 'ISO'/'SAE': use ISO or SAE orientation of the normal for offset computation
    %>
    %> **Output**
    %>
    %> - `XY`: matrix `2 x n` of the evaluated points
    %> - `X`: vector of the x-coordinates of the evaluated point derivatives
    %> - `Y`: vector of the y-coordinates of the evaluated point derivatives
    %>
    function varargout = eval_DDD( self, varargin )
      [ varargout{1:nargout} ] = ...
        feval( self.mexName, 'eval_DDD', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the angle of the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   theta = ref.theta( s );
    %>
    %> ```
    %>
    function th = theta(self, s)
      th = feval( self.mexName, 'theta', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the angle derivatives (curvature) of the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   theta = ref.theta_D( s );
    %>
    %> ```
    %>
    function th = theta_D(self, s)
      % evaluate angle derivative [curvature] at `s`
      th = feval( self.mexName, 'theta_D', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the angle second derivatuve of the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   theta = ref.theta_DD( s );
    %>
    %> ```
    %>
    function th = theta_DD(self, s)
      % evaluate angle second derivative at `s`
      th = feval( self.mexName, 'theta_DD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the angle third derivative of the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   theta = ref.theta_DDD( s );
    %>
    %> ```
    %>
    function th = theta_DDD(self, s)
      % evaluate angle third derivative at `s`
      th = feval( self.mexName, 'theta_DDD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the curvature of the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   theta = ref.kappa( s );
    %>
    %> ```
    %>
    function th = kappa(self, s)
      % evaluate curvature at `s`
      th = feval( self.mexName, 'kappa', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the curvature derivative of the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   theta = ref.kappa_D( s );
    %>
    %> ```
    %>
    function th = kappa_D(self, s)
      % evaluate curvature derivative at `s`
      th = feval( self.mexName, 'kappa_D', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the curvature second derivative of the curve at curvilinear coordinate `s`.
    %> Argument `s` may be a vector for multiple evaluations.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   theta = ref.kappa_DD( s );
    %>
    %> ```
    %>
    %>
    function th = kappa_DD(self, s)
      % evaluate curvature second derivative at `s`
      th = feval( self.mexName, 'kappa_DD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate initial point of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   [ x0, y0 ] = ref.xy_begin();
    %>
    %> ```
    %>
    function [ x, y ] = xy_begin( self )
      [ x, y ] = feval( self.mexName, 'xy_begin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function [ x, y ] = xyBegin( self )
      [ x, y ] = self.xy_begin();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate final point of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   [ x1, y1 ] = ref.xy_end();
    %>
    %> ```
    %>
    function [ x, y ] = xy_end( self )
      [ x, y ] = feval( self.mexName, 'xy_end', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function [ x, y ] = xyEnd( self )
      [ x, y ] = self.xy_end();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate initial x-coordinate of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   x0 = ref.x_begin();
    %>
    %> ```
    %>
    function X0 = x_begin( self )
      X0 = feval( self.mexName, 'x_begin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function X0 = xBegin( self )
      X0 = self.x_begin();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate final x-coordinate of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   x1 = ref.x_end();
    %>
    %> ```
    %>
    function X1 = x_end( self )
      X1 = feval( self.mexName, 'x_end', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function X1 = xEnd( self )
      X1 = self.x_end();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate initial y-coordinate of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   y0 = ref.y_begin();
    %>
    %> ```
    %>
    function Y0 = y_begin( self )
      Y0 = feval( self.mexName, 'y_begin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function Y0 = yBegin( self )
      Y0 = self.y_begin();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate final y-coordinate of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   y1 = ref.y_end();
    %>
    %> ```
    %>
    function Y1 = y_end( self )
      Y1 = feval( self.mexName, 'y_end', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function Y1 = yEnd( self )
      Y1 = self.y_end();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate initial angle of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   theta = ref.theta_begin();
    %>
    %> ```
    %>
    function th0 = theta_begin( self )
      th0 = feval( self.mexName, 'theta_begin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function th0 = thetaBegin( self )
      th0 = self.theta_begin();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate final angle of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   theta = ref.theta_end();
    %>
    %> ```
    %>
    function th1 = theta_end( self )
      th1 = feval( self.mexName, 'theta_end', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function th1 = thetaEnd( self )
      th1 = self.theta_end();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate initial curvature of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   kappa0 = ref.kappa_begin();
    %>
    %> ```
    %>
    function kappa0 = kappa_begin( self )
      kappa0 = feval( self.mexName, 'kappa_begin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function kappa0 = kappaBegin( self )
      kappa0 = self.kappa_begin();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate final curvature of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   kappa1 = ref.kappa_end();
    %>
    %> ```
    %>
    function kappa1 = kappa_end( self )
      kappa1 = feval( self.mexName, 'kappa_end', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function kappa1 = kappaEnd( self )
      kappa1 = self.kappa_end();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Return the length of the curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   length = ref.length();
    %>
    %> ```
    %>
    function res = length( self, varargin )
      res = feval( self.mexName, 'length', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ p1, p2 ] = points( self )
      [ p1, p2 ] = feval( self.mexName, 'points', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the bounding box triangles of curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   [P1,P2,P3] = ref.bbTriangles();
    %>   [P1,P2,P3] = ref.bbTriangles(max_angle,max_size);
    %>   [P1,P2,P3] = ref.bbTriangles(max_angle,max_size,offs);
    %>   [P1,P2,P3] = ref.bbTriangles(max_angle,max_size,offs,'ISO');
    %>   [P1,P2,P3] = ref.bbTriangles(max_angle,max_size,offs,'SAE');
    %>
    %> ```
    %>
    %> **Optional Arguments**
    %>
    %> - `max_angle`: maximum curve angle variation admitted in a triangle
    %> - `max_size`: maximum triangles size
    %> - `offs`: offset of the curve used in computation
    %> - 'ISO'/'SAE': use ISO or SAE orientation of the normal for the offset
    %>
    %> **Output**
    %>
    %> - `P1`: `2 x n` matrix with the first points of the triangles
    %> - `P2`: `2 x n` matrix with the second points of the triangles
    %> - `P3`: `2 x n` matrix with the third points of the triangles
    %>
    function [P1,P2,P3] = bbTriangles( self, varargin )
      [P1,P2,P3] = feval( self.mexName, 'bbTriangles', ...
        self.objectHandle, varargin{:} ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the point at minimum distance of another point on the curve.
    %> `qx` and `qy` may be vectors so that the return values are vectors too.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   [ x, y, s, t, iflag, dst ] = ref.closest_point( qx, qy );
    %>   [ x, y, s, t, iflag, dst ] = ref.closest_point( qx, qy, offs );
    %>   [ x, y, s, t, iflag, dst ] = ref.closest_point( qx, qy, offs, 'ISO' );
    %>   [ x, y, s, t, iflag, dst ] = ref.closest_point( qx, qy, offs, 'SAE' );
    %>
    %> ```
    %>
    %> **Optional Arguments**
    %>
    %> - offs: offset of the curve used in computation
    %> - 'ISO'/'SAE': use ISO or SAE orientation of the normal for the offset
    %>
    %> **Output**
    %>
    %> - `x`, `y`: Point at minimum distance from `(qx,qy)` on the curve.
    %> - `s`, `t`: Curvilinear coordinates of the point `(qx,qy)`.
    %> - `iflag`: `iflag < 0` some error in computation, iflag >0 is the numer of segment
    %>   containing the point at minimum distance.
    %> - `dst`: point curve distance.
    %>
    function varargout = closest_point( self, qx, qy, varargin )
      % iflag > 0 = segment number
      [ varargout{1:nargout} ] = feval( self.mexName, 'closest_point', ...
        self.objectHandle, qx, qy, varargin{:} ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function varargout = closestPoint( self, qx, qy, varargin )
      % iflag > 0 = segment number
      [ varargout{1:nargout} ] = self.closest_point( qx, qy, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Evaluate the distance of a point `(qx,qy)` to the curve.
    %> `qx` and `qy` may be vectors so that the return values are vectors too.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   dst = ref.distance( qx, qy );
    %>   dst = ref.distance( qx, qy, offs );
    %>   dst = ref.distance( qx, qy, offs, 'ISO' );
    %>   dst = ref.distance( qx, qy, offs, 'SAE' );
    %>
    %> ```
    %>
    %> **Optional Arguments**
    %>
    %> - `offs`: offset of the curve used in computation
    %> - 'ISO'/'SAE': use ISO or SAE orientation of the normal for the offset
    %>
    function dst = distance( self, qx, qy, varargin )
      dst = feval( self.mexName, 'distance', ...
        self.objectHandle, qx, qy, varargin{:} ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Check if two curve collide.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   ok = ref.collision( obj );
    %>   ok = ref.collision( obj, offs, offs1 );
    %>   ok = ref.collision( obj, offs, offs1, 'ISO' );
    %>   ok = ref.collision( obj, offs, offs1, 'SAE' );
    %>
    %> ```
    %>
    %> **Optional Arguments**
    %>
    %> - `offs`, `offs1`: offset of the curves used in computation
    %> - 'ISO'/'SAE': use ISO or SAE orientation of the normal for the offsets
    %>
    function ok = collision( self, OBJ, varargin )
      ok = feval( self.mexName, 'collision', ...
        self.objectHandle, OBJ.obj_handle(), OBJ.is_type(), ...
        varargin{:} ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Intersect two curves.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   [s1,s2] = ref.intersect( obj );
    %>   [s1,s2] = ref.intersect( obj, offs, offs1 );
    %>   [s1,s2] = ref.intersect( obj, offs, offs1, 'ISO' );
    %>   [s1,s2] = ref.intersect( obj, offs, offs1, 'SAE' );
    %>
    %> ```
    %>
    %> - `s1`: curvilinear coordinates of the intersections on the first curve
    %> - `s2`: curvilinear coordinates of the intersections on the second curve
    %>
    %> **Optional Argument**
    %>
    %> - `offs`, `offs1`: offset of the curves used in computation
    %> - 'ISO'/'SAE': use ISO or SAE orientation of the normal for the offsets
    %>
    function [s1,s2] = intersect( self, OBJ, varargin )
      [s1,s2] = feval( self.mexName, 'intersect', ...
        self.objectHandle, OBJ.obj_handle(), OBJ.is_type(), ...
        varargin{:} ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Print on the console some information on the stored curve.
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>
    %>   ref.info();
    %>
    %> ```
    %>
    function info( self )
      feval( self.mexName, 'info', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get the curvilinear coordinates of the point `(x,y)`
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   [s,t] = ref.find_coord( x, y );
    %>
    %> ```
    %>
    %> - `s`: curvilinear coordinate along the curve
    %> - `t`: curvilinear coordinate along the normal of the curve
    %>
    function [s,t] = find_coord( self, x, y )
      [s,t] = feval( self.mexName, 'findST', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Activate the use of AABB three in intersection/collision computations
    %>
    function yesAABBtree( self )
      feval( self.mexName, 'yesAABBtree' );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Deactivate the use of AABB three in intersection/collision computations
    %>
    function noAABBtree( self )
      feval( self.mexName, 'noAABBtree' );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot a triangle BBOX
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plot_tbox( P1, P2, P3 );
    %>   ref.plot_tbox( P1, P2, P3, 'Color', 'red' );
    %>
    %> ```
    %>
    %>
    function plot_tbox( self, P1, P2, P3, varargin )
      for k=1:size(P1,2)
        pp1 = P1(:,k);
        pp2 = P2(:,k);
        pp3 = P3(:,k);
        plot( [pp1(1),pp2(1),pp3(1),pp1(1)], ...
              [pp1(2),pp2(2),pp3(2),pp1(2)], varargin{:} );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function plotTBox( self, P1, P2, P3, varargin )
      self.plot_tbox( P1, P2, P3, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the bounding box of the curve
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plot_bbox();
    %>   ref.plot_bbox( 'Color', 'red' );
    %>
    %> ```
    %>
    function plot_bbox( self, varargin )
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
    %>
    %> \deprecated  whill be removed in future version
    %>
    function plotBBox( self, varargin )
      self.plot_bbox( varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the covering triangles of the curve
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plot_triangles()
    %>   ref.plot_triangles('red','FaceAlpha', 0.5);
    %>
    %> ```
    %>
    function plot_triangles( self, varargin )
      [p1,p2,p3] = self.bbTriangles();
      for k=1:size(p1,2)
        x = [ p1(1,k), p2(1,k), p3(1,k), p1(1,k) ];
        y = [ p1(2,k), p2(2,k), p3(2,k), p1(2,k) ];
        %fill( x, y, 'red','FaceAlpha', 0.5 );
        fill( x, y, varargin{:});
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function plotTriangles( self, varargin )
      self.plot_triangles( varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
