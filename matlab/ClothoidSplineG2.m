classdef ClothoidSplineG2 < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    x;
    y;
    theta;
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
      self.theta     = (phi_L./len_L+phi_R./len_R)./(1./len_L+1./len_R) ;
      self.theta_min = max(phi_L,phi_R)-0.99*pi ;
      self.theta_max = min(phi_L,phi_R)+0.99*pi ;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function k = KAPPA( ~, theta0, theta )
      x = theta0.^2 ;
      a = -3.714 + x * 0.178 ;
      b = -1.913 - x * 0.0753 ;
      c =  0.999 + x * 0.03475 ;
      d =  0.191 - x * 0.00703 ;
      e =  0.500 - x * -0.00172 ;
      k = a.*theta0+b.*theta+c.*(d.*theta0+e.*theta).^3 ;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
    function f = TG_objective( ~, ~ )
      f = 0 ;
    end
    %
    function g = TG_gradient( ~, theta )
      N = length(theta) ;
      g = zeros(N,1) ;
    end
    %
    function F = TG_constraints( self, theta )
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
      F  = [ kL(1:end-1)-k(2:end) ; ...
             self.diff2pi( theta(1)   - self.theta_I ) ; ...
             self.diff2pi( theta(end) - self.theta_F ) ] ;
    end
    % 
    function JAC = TG_jacobian( self, theta )
      N  = length(theta) ;
      L  = zeros(1,N-1) ; L_1  = zeros(1,N-1) ; L_2  = zeros(1,N-1) ;
      k  = zeros(1,N-1) ; k_1  = zeros(1,N-1) ; k_2  = zeros(1,N-1) ;
      dk = zeros(1,N-1) ; dk_1 = zeros(1,N-1) ; dk_2 = zeros(1,N-1) ;
      c  = ClothoidCurve() ;

      JAC = sparse(N,N) ;
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
      for j=1:N-2
        JAC(j,j)   =  k_1(j) + dk_1(j)*L(j) + dk(j)*L_1(j) ;
        JAC(j,j+1) =  k_2(j) + dk_2(j)*L(j) + dk(j)*L_2(j) - k_1(j+1) ;
        JAC(j,j+2) = -k_2(j+1) ;
      end
      JAC(N-1,1) = 1 ;
      JAC(N,N)   = 1 ;
    end
    %
    function JAC = TG_jacobianstructure( self )
      N  = length(self.theta) ;
      JAC = sparse(N,N) ;
      for j=1:N-2
        JAC(j,j)   = 1 ;
        JAC(j,j+1) = 1 ;
        JAC(j,j+2) = 1 ;
      end
      JAC(N-1,1)   = 1 ;
      JAC(N,N)     = 1 ;
    end
  end

  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function self = ClothoidSplineG2()
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function delete( self )
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function clots = buildTG( self, x, y )
      %
      % copy data
      %
      self.x = x ;
      self.y = y ;
      %
      % Compute guess angles
      %
      self.guess_angle()
      self.theta_I = self.theta(1) ;
      self.theta_F = self.theta(end) ;
      
      %
      % Compute solution using IPOPT
      %
      N = length(self.theta) ;
      options = {} ;

      options.ub = self.theta_max ;
      options.lb = self.theta_min ;

      % The constraint functions are bounded to zero
      options.cl = zeros(N,1) ; %  constraints
      options.cu = zeros(N,1) ;

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
      funcs.objective         = @(theta) self.TG_objective(theta);
      funcs.constraints       = @(theta) self.TG_constraints(theta);
      funcs.gradient          = @(theta) self.TG_gradient(theta);
      funcs.jacobian          = @(theta) self.TG_jacobian(theta);
      funcs.jacobianstructure = @() self.TG_jacobianstructure();
      
      %options.ipopt.jacobian_approximation = 'finite-difference-values';    
      options.ipopt.hessian_approximation  = 'limited-memory';
      %options.ipopt.limited_memory_update_type = 'bfgs' ; % {bfgs}, sr1 = 6; % {6}
      %options.ipopt.limited_memory_update_type = 'sr1' ;
      options.ipopt.limited_memory_update_type = 'bfgs' ; % {bfgs}, sr1 = 6; % {6}
  
      tic
      [theta, info] = ipopt(self.theta,funcs,options);
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
        clots.push_back( x(j),   y(j),   theta(j), ...
                         x(j+1), y(j+1), theta(j+1) ) ;
      end
      clots.plot(0.01);
      %stats.nevalF    = output.funcCount ;  
      %stats.iter      = output.iterations ;
      %stats.Fvalue    = sqrt(resnorm/N) ;
      %stats.Fgradnorm = output.firstorderopt ;
    end
    % 
  end
end