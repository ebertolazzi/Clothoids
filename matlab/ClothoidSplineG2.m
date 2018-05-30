classdef ClothoidSplineG2 < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    t_type;
    x;
    y;
    theta_guess;
    theta_I;
    theta_F;
    theta_min;
    theta_max;
  end

  methods (Hidden = true)

    function guess_angle( self )
      %
      % Compute guess angles
      %
      NPTS   = length(self.x) ;
      phi    = zeros(1,NPTS-1) ;
      len    = zeros(1,NPTS-1) ;
      d      = [ self.x(2)-self.x(1) ; self.y(2)-self.y(1) ] ;
      phi(1) = atan2(d(2),d(1)) ;
      len(1) = norm(d,2) ;
      for k=2:NPTS-1
        d      = [ self.x(k+1)-self.x(k) ; self.y(k+1)-self.y(k) ] ;
        len(k) = norm(d,2) ;
        phi(k) = atan2(d(2),d(1)) ;
        df     = phi(k)-phi(k-1) ;
        df     = df - round(df/(2*pi))*2*pi ;
        phi(k) = phi(k-1)+df ;
      end
      phi_L = [ phi(1) phi ] ; phi_R = [ phi phi(end) ] ;
      len_L = [ len(1) len ] ; len_R = [ len len(end) ] ;
      %theta     = (phi_L+phi_R)./2 ;
      %theta     = (phi_L.*len_L+phi_R.*len_R)./(len_L+len_R) ;
      self.theta_guess = (phi_L./len_L+phi_R./len_R)./(1./len_L+1./len_R) ;
      self.theta_min   = max(phi_L,phi_R)-0.99*pi ;
      self.theta_max   = min(phi_L,phi_R)+0.99*pi ;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function res = diff2pi( ~, in )
      res = in-2*pi*round(in/(2*pi)) ;
    end
    %  ___ ____   ___  ____ _____
    % |_ _|  _ \ / _ \|  _ \_   _|
    %  | || |_) | | | | |_) || |
    %  | ||  __/| |_| |  __/ | |
    % |___|_|    \___/|_|    |_|
    %
    function f = objective( self, theta )
      cL = ClothoidCurve() ;
      cR = ClothoidCurve() ;
      cL.build_G1( self.x(1), self.y(1), theta(1), ...
                   self.x(2), self.y(2), theta(2) ) ;
      cR.build_G1( self.x(end-1), self.y(end-1), theta(end-1), ...
                   self.x(end),   self.y(end),   theta(end) ) ;
      N = length(theta) ;

      switch self.t_type
      case 'P1'
        f = 0 ;
      case 'P2'
        %[kL,~,~] = clothoid(1,theta) ;
        %[kR,dkR,LR] = clothoid(N-1,theta) ;
        f = self.diff2pi(theta(1)-theta(end))^2; % (kR+LR*dkR-kL)^2 ;
      case 'P4'
        f = cL.kappa_D()^2+cR.kappa_D()^2 ;
      case 'P5'
        f = cL.length()+cR.length();
      case 'P6'
        c = ClothoidCurve() ;
        f = 0 ;
        for j=1:N-1
          c.build_G1( self.x(j),   self.y(j),   theta(j), ...
                      self.x(j+1), self.y(j+1), theta(j+1) ) ;
          f = f + c.length();
        end
      case 'P7'
        c = ClothoidCurve() ;
        f = 0 ;
        for j=1:N-1
          c.build_G1( self.x(j),   self.y(j),   theta(j), ...
                      self.x(j+1), self.y(j+1), theta(j+1) ) ;
          [~,~,~,k,dk,L] = c.getPars();
          f = f + (dk^2*L^3)/3 + dk*L^2*k + L*k^2 ;
        end
      case 'P8'
        c = ClothoidCurve() ;
        f = 0 ;
        for j=1:N-1
          c.build_G1( self.x(j),   self.y(j),   theta(j), ...
                      self.x(j+1), self.y(j+1), theta(j+1) ) ;
          [~,~,~,~,dk,L] = c.getPars();
          f = f + L*dk^2 ;
        end
      case 'P9'
        c = ClothoidCurve() ;
        f = 0 ;
        for j=1:N-1
          c.build_G1( self.x(j),   self.y(j),   theta(j), ...
                      self.x(j+1), self.y(j+1), theta(j+1) ) ;
          [~,~,~,k,dk,L] = c.getPars();
          f = f + (dk^4*L^5)/5+k*dk^3*L^4+(2*L^3*k^2+L)*dk^2+2*k^3*dk*L^2+k^4*L ;
        end
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = gradient( self, theta )
      cL = ClothoidCurve() ;
      cR = ClothoidCurve() ;
      [LL_D, ~, dkL_D] = cL.build_G1( self.x(1), self.y(1), theta(1), ...
                                      self.x(2), self.y(2), theta(2) ) ;
      [LR_D, ~, dkR_D] = cR.build_G1( self.x(end-1), self.y(end-1), theta(end-1), ...
                                      self.x(end),   self.y(end),   theta(end) ) ;

      dkL = cL.kappa_D() ;
      dkR = cR.kappa_D() ;

      N = length(theta) ;
      g = zeros(N,1) ;
      switch self.t_type
      case 'P2'
        g(1)   = + 2*self.diff2pi(theta(1)-theta(end)) ;
        g(end) = - 2*self.diff2pi(theta(1)-theta(end)) ;
      case 'P4'
        g(1)     = 2*dkL*dkL_D(1) ;
        g(2)     = 2*dkL*dkL_D(2) ;
        g(end-1) = 2*dkR*dkR_D(1) ;
        g(end)   = 2*dkR*dkR_D(2) ;
      case 'P5'
        g(1)     = LL_D(1) ;
        g(2)     = LL_D(2) ;
        g(end-1) = LR_D(1) ;
        g(end)   = LR_D(2) ;
      case 'P6'
        c = ClothoidCurve() ;
        for j=1:N-1
          [L_D,~,~] = c.build_G1( self.x(j),   self.y(j),   theta(j), ...
                                  self.x(j+1), self.y(j+1), theta(j+1) ) ;
          g(j)   = g(j)   + L_D(1) ;
          g(j+1) = g(j+1) + L_D(2) ;
        end
      case 'P7'
        c = ClothoidCurve() ;
        for j=1:N-1
          [L_D, k_D, dk_D] = c.build_G1( self.x(j),   self.y(j),   theta(j), ...
                                         self.x(j+1), self.y(j+1), theta(j+1) ) ;
          [~,~,~,k,dk,L] = c.getPars();
          g(j)   = g(j)   + 2*(dk*dk_D(1)*L^3)/3 + (dk^2*L^2*L_D(1)) ...
                          + dk_D(1)*L^2*k + dk*2*L*L_D(1)*k + dk*L^2*k_D(1) ...
                          + L_D(1)*k^2 + 2*L*k*k_D(1) ;
          g(j+1) = g(j+1) + 2*(dk*dk_D(2)*L^3)/3 + (dk^2*L^2*L_D(2)) ...
                          + dk_D(2)*L^2*k + dk*2*L*L_D(2)*k + dk*L^2*k_D(2) ...
                          + L_D(2)*k^2 + 2*L*k*k_D(2) ;
        end
      case 'P8'
        c = ClothoidCurve() ;
        for j=1:N-1
          [L_D,~,dk_D] = c.build_G1( self.x(j),   self.y(j),   theta(j), ...
                                         self.x(j+1), self.y(j+1), theta(j+1) ) ;
          [~,~,~,~,dk,L] = c.getPars();
          g(j)   = g(j)   + 2*L*dk*dk_D(1) + L_D(1)*dk^2 ;
          g(j+1) = g(j+1) + 2*L*dk*dk_D(2) + L_D(2)*dk^2 ;
        end
      case 'P9'
        c = ClothoidCurve() ;
        for j=1:N-1
          [L_D, k_D, dk_D] = c.build_G1( self.x(j),   self.y(j),   theta(j), ...
                                         self.x(j+1), self.y(j+1), theta(j+1) ) ;
          [~,~,~,k,dk,L] = c.getPars();
          A = dk^4*L^4+4*k*dk^3*L^3+(6*L^2*k^2+1)*dk^2+4*k^3*dk*L+k^4 ;
          B = (4*L^3*k^2+2*L)*dk+3*k*dk^2*L^4+(4/5)*dk^3*L^5+2*k^3*L^2 ; 
          C = L^4*dk^3+4*L^3*dk^2*k+6*L^2*dk*k^2+4*L*k^3 ;
          g(j)   = g(j)   + A*L_D(1) + B*dk_D(1) + C*k_D(1) ;
          g(j+1) = g(j+1) + A*L_D(2) + B*dk_D(2) + C*k_D(2) ;
        end
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = constraints( self, theta )
      N  = length(theta) ;
      k  = zeros(N-1,1) ;
      dk = zeros(N-1,1) ;
      L  = zeros(N-1,1) ;
      c  = ClothoidCurve() ;
      for j=1:N-1
        c.build_G1( self.x(j),   self.y(j),   theta(j), ...
                    self.x(j+1), self.y(j+1), theta(j+1) ) ;
        k(j)  = c.kappaBegin() ;
        dk(j) = c.kappa_D() ;
        L(j)  = c.length() ;
      end
      kL = k+dk.*L ;
      switch self.t_type
      case 'P1'
        F  = [ kL(1:end-1)-k(2:end) ; ...
               self.diff2pi( theta(1)   - self.theta_I ) ; ...
               self.diff2pi( theta(end) - self.theta_F ) ] ;
      case 'P2'
        F = [ kL(1:end-1)-k(2:end) ; kL(end)-k(1) ] ;
      otherwise
        F = kL(1:end-1)-k(2:end) ;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function JAC = jacobian( self, theta )
      N  = length(theta) ;
      L  = zeros(1,N-1) ; L_1  = zeros(1,N-1) ; L_2  = zeros(1,N-1) ;
      k  = zeros(1,N-1) ; k_1  = zeros(1,N-1) ; k_2  = zeros(1,N-1) ;
      dk = zeros(1,N-1) ; dk_1 = zeros(1,N-1) ; dk_2 = zeros(1,N-1) ;
      c  = ClothoidCurve() ;

      for j=1:N-1
        [L_D,k_D,dk_D] = c.build_G1( self.x(j),   self.y(j),   theta(j), ...
                                     self.x(j+1), self.y(j+1), theta(j+1) ) ;
        k(j)    = c.kappaBegin() ;
        dk(j)   = c.kappa_D() ;
        L(j)    = c.length() ;
        L_1(j)  = L_D(1)  ; L_2(j)  = L_D(2) ; 
        k_1(j)  = k_D(1)  ; k_2(j)  = k_D(2) ; 
        dk_1(j) = dk_D(1) ; dk_2(j) = dk_D(2) ; 
      end

      switch self.t_type
      case 'P1'
        JAC = sparse(N,N) ;
      case 'P2'
        JAC = sparse(N-1,N) ;
      otherwise
        JAC = sparse(N-2,N) ;
      end

      for j=1:N-2
        JAC(j,j)   =  k_1(j) + dk_1(j)*L(j) + dk(j)*L_1(j) ;
        JAC(j,j+1) =  k_2(j) + dk_2(j)*L(j) + dk(j)*L_2(j) - k_1(j+1) ;
        JAC(j,j+2) = -k_2(j+1) ;
      end

      switch self.t_type
      case 'P1'
        JAC(N-1,1) = 1 ;
        JAC(N,N)   = 1 ;
      case 'P2'
        JAC(N-1,1)   = -k_1(1) ;
        JAC(N-1,2)   = -k_2(1) ;
        JAC(N-1,N-1) = k_1(end)+L_1(end)*dk(end)+L(end)*dk_1(end) ;
        JAC(N-1,N)   = k_2(end)+L_2(end)*dk(end)+L(end)*dk_2(end) ;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function JAC = jacobianstructure( self )
      N = length(self.theta_guess) ;
      
      switch self.t_type
      case 'P1'
        JAC = sparse(N,N) ;
      case 'P2'
        JAC = sparse(N-1,N) ;
      otherwise
        JAC = sparse(N-2,N) ;
      end
      
      for j=1:N-2
        JAC(j,j)   = 1 ;
        JAC(j,j+1) = 1 ;
        JAC(j,j+2) = 1 ;
      end
      
      switch self.t_type
      case 'P1'
        JAC(N-1,1) = 1 ;
        JAC(N,N)   = 1 ;
      case 'P2'
        JAC(N-1,1)   = 1;
        JAC(N-1,2)   = 1;
        JAC(N-1,N-1) = 1;
        JAC(N-1,N)   = 1;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function clots = build_internal( self, x, y, varargin )
      %
      % copy data
      %
      self.x = x ;
      self.y = y ;
      %
      % Compute guess angles
      %
      self.guess_angle()
      if nargin == 5
        self.theta_I = varargin{1} ;
        self.theta_F = varargin{2} ;
      else
        self.theta_I = self.theta_guess(1) ;
        self.theta_F = self.theta_guess(end) ;
      end
      %
      % Compute solution using IPOPT
      %
      N = length(self.theta_guess) ;
      options = {} ;

      options.ub = self.theta_max ;
      options.lb = self.theta_min ;
                  
      switch self.t_type
      case 'P1'
        nc = N ;
      case 'P2'
        nc = N-1 ;
      otherwise
        nc = N-2 ;
      end

      % The constraint functions are bounded to zero
      options.cl = zeros(nc,1) ; %  constraints
      options.cu = zeros(nc,1) ;

      % Set the IPOPT options.
      options.ipopt.jac_d_constant      = 'no' ;
      options.ipopt.hessian_constant    = 'no' ;
      options.ipopt.mu_strategy         = 'adaptive';
      options.ipopt.max_iter            = 400 ;
      options.ipopt.tol                 = 1e-10 ;
      options.ipopt.derivative_test_tol = 1e-4 ;
      options.ipopt.derivative_test     = 'first-order' ; %% default 'none'
      options.ipopt.derivative_test_perturbation = 1e-8 ;

      % The callback functions.
      funcs.objective         = @(theta) self.objective(theta);
      funcs.constraints       = @(theta) self.constraints(theta);
      funcs.gradient          = @(theta) self.gradient(theta);
      funcs.jacobian          = @(theta) self.jacobian(theta);
      funcs.jacobianstructure = @()      self.jacobianstructure();
      
      %options.ipopt.jacobian_approximation = 'finite-difference-values';    
      options.ipopt.hessian_approximation  = 'limited-memory';
      %options.ipopt.limited_memory_update_type = 'bfgs' ; % {bfgs}, sr1 = 6; % {6}
      %options.ipopt.limited_memory_update_type = 'sr1' ;
      options.ipopt.limited_memory_update_type = 'bfgs' ; % {bfgs}, sr1 = 6; % {6}
  
      tic
      [theta, info] = ipopt(self.theta_guess,funcs,options);
      stats.elapsed = toc ;

      info;

      %options = optimset(varargin{:});
      %[theta,resnorm,~,~,output,~,~] = lsqnonlin( @target, theta, [], [], options ) ;
      %
      % Compute spline parameters
      %
      clots = ClothoidList() ;
      N     = length(theta) ;
      clots.reserve(N-1);
      for j=1:N-1
        clots.push_back( x(j), y(j), theta(j), x(j+1), y(j+1), theta(j+1) ) ;
      end
    end
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = ClothoidSplineG2()
      self.t_type = 'P1' ;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( ~ )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function clots = buildP1( self, x, y )
      self.t_type = 'P1';
      clots = self.build_internal( x, y ) ;
    end 
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function clots = buildP2( self, x, y )
      self.t_type = 'P2';
      clots = self.build_internal( x, y ) ;
    end 
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function clots = buildP4( self, x, y )
      self.t_type = 'P4';
      clots = self.build_internal( x, y ) ;
    end 
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function clots = buildP5( self, x, y )
      self.t_type = 'P5';
      clots = self.build_internal( x, y ) ;
    end 
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function clots = buildP6( self, x, y )
      self.t_type = 'P6';
      clots = self.build_internal( x, y ) ;
    end 
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function clots = buildP7( self, x, y )
      self.t_type = 'P7';
      clots = self.build_internal( x, y ) ;
    end 
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function clots = buildP8( self, x, y )
      self.t_type = 'P8';
      clots = self.build_internal( x, y ) ;
    end 
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function clots = buildP9( self, x, y )
      self.t_type = 'P9';
      clots = self.build_internal( x, y ) ;
    end 
    % 
  end
end