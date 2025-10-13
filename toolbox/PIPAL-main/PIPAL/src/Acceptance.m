% Author      : Frank E. Curtis
% Description : Class for line search quantities and routines.

% Acceptance class
classdef Acceptance < handle
  
  % Class properties (private set access)
  properties (SetAccess = private)
    
    p0 = 0; % Fraction-to-the-boundary steplength
    p  = 0; % Primal steplength
    d  = 0; % Dual steplength
    s  = 0; % Bool for second-order correction
    
  end
  
  % Class methods
  methods
    
    % Backtracking line search
    function backtracking(a,p,i,c,z,d)
      
      % Store current values
      x = z.x; f = z.f;
      cE = z.cE; r1 = z.r1; r2 = z.r2; lE = z.lE;
      cI = z.cI; s1 = z.s1; s2 = z.s2; lI = z.lI;
      phi = z.phi;
      
      % Backtracking loop
      while a.p >= eps
        
        % Set trial point
        z.updatePoint  (i,  d,a);
        z.evalFunctions(i,c    );
        
        % Check for function evaluation error
        if z.err == 0
          
          % Set remaining trial values
          z.evalSlacks(p,i);
          z.evalMerit (  i);
          
          % Check for nonlinear fraction-to-boundary violation
          ftb = 0;
          if i.nE > 0, ftb = ftb + sum(z.r1<min(p.ls_frac,z.mu)*r1) + sum(z.r2<min(p.ls_frac,z.mu)*r2); end;
          if i.nI > 0, ftb = ftb + sum(z.s1<min(p.ls_frac,z.mu)*s1) + sum(z.s2<min(p.ls_frac,z.mu)*s2); end;
          
          % Check Armijo condition
          if ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*max(d.qtred,0)
            
            % Reset variables and return
            z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,z.f,z.cE,z.cI,z.phi); return;
            
          else
            
            % Reset variables
            z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,f,cE,cI,phi);
            
          end;
          
        else
          
          % Reset variables
          z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,f,cE,cI,phi);
          
        end
        
        % Reduce steplength
        a.p = p.ls_factor*a.p;
        
      end
      
    end
    
    % Fraction-to-boundary line search
    function fractionToBoundary(a,p,i,z,d)
      
      % Initialize primal fraction-to-boundary
      a.p0 = 1;
      
      % Update primal fraction-to-boundary for constraint slacks
      if i.nE > 0, a.p0 = min([a.p0; (min(p.ls_frac,z.mu)-1)*z.r1(d.r1<0)./d.r1(d.r1<0); (min(p.ls_frac,z.mu)-1)*z.r2(d.r2<0)./d.r2(d.r2<0)]); end;
      if i.nI > 0, a.p0 = min([a.p0; (min(p.ls_frac,z.mu)-1)*z.s1(d.s1<0)./d.s1(d.s1<0); (min(p.ls_frac,z.mu)-1)*z.s2(d.s2<0)./d.s2(d.s2<0)]); end;
      
      % Initialize primal steplength
      a.p = a.p0;
      
      % Initialize dual fraction-to-boundary
      a.d = 1;
      
      % Update dual fraction-to-boundary for constraint multipliers
      if i.nE > 0, a.d = min([a.d; (min(p.ls_frac,z.mu)-1)*(1+z.lE(d.lE<0))./d.lE(d.lE<0); (1-min(p.ls_frac,z.mu))*(1-z.lE(d.lE>0))./d.lE(d.lE>0)]); end;
      if i.nI > 0, a.d = min([a.d; (min(p.ls_frac,z.mu)-1)*(0+z.lI(d.lI<0))./d.lI(d.lI<0); (1-min(p.ls_frac,z.mu))*(1-z.lI(d.lI>0))./d.lI(d.lI>0)]); end;
      
    end
    
    % Full step search for trial penalty parameters
    function b = fullStepCheck(a,p,i,c,z,d)
      
      % Initialize boolean
      b = 0;
      
      % Set current and last penalty parameters
      rho      = z.rho;
      rho_temp = z.rho_;
      
      % Loop through last penalty parameters
      while rho < rho_temp
        
        % Set penalty parameter
        z.setRho(rho_temp);
        
        % Evaluate merit
        z.evalMerit(i);
        
        % Store current values
        x = z.x; f = z.f;
        cE = z.cE; r1 = z.r1; r2 = z.r2; lE = z.lE;
        cI = z.cI; s1 = z.s1; s2 = z.s2; lI = z.lI;
        phi = z.phi;
        
        % Set trial point
        z.updatePoint  (i,  d,a);
        z.evalFunctions(i,c    );
        
        % Check for function evaluation error
        if z.err == 0
          
          % Set remaining trial values
          z.evalSlacks(p,i);
          z.evalMerit (  i);
          
          % Check for nonlinear fraction-to-boundary violation
          ftb = 0;
          if i.nE > 0, ftb = ftb + sum(z.r1<min(p.ls_frac,z.mu)*r1) + sum(z.r2<min(p.ls_frac,z.mu)*r2); end;
          if i.nI > 0, ftb = ftb + sum(z.s1<min(p.ls_frac,z.mu)*s1) + sum(z.s2<min(p.ls_frac,z.mu)*s2); end;
          
          % Check Armijo condition
          if ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*max(d.qtred,0)
            
            % Reset variables, set boolean, and return
            z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,z.f,z.cE,z.cI,z.phi); b = 1; return;
            
          else
            
            % Reset variables
            z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,f,cE,cI,phi);
            
          end
          
        else
          
          % Reset variables
          z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,f,cE,cI,phi);
          
        end
        
        % Decrease rho
        rho_temp = p.rho_factor*rho_temp;
        
      end
      
      % Set rho
      z.setRho(rho);
      
      % Evaluate merit
      z.evalMerit(i);
      
    end
    
    % Line search
    function lineSearch(a,p,i,c,z,d)
      
      % Check fraction-to-boundary rule
      a.fractionToBoundary(p,i,z,d);
      
      % Check for full step for trial penalty parameters
      b = a.fullStepCheck(p,i,c,z,d);
      
      % Run second-order correction
      a.s = 0; if b == 0, b = a.secondOrderCorrection(p,i,c,z,d); if b == 2, a.s = 1; end; end;
      
      % Run backtracking line search
      if b == 0, a.backtracking(p,i,c,z,d); end;
      
    end
    
    % Second-order Correction
    function b = secondOrderCorrection(a,p,i,c,z,d)
      
      % Initialize flag
      b = 0;
      
      % Store current iterate values
      x = z.x; f = z.f;
      cE = z.cE; r1 = z.r1; r2 = z.r2; lE = z.lE;
      cI = z.cI; s1 = z.s1; s2 = z.s2; lI = z.lI;
      phi = z.phi; v = z.v;
      
      % Set trial point
      z.updatePoint  (i,  d,a);
      z.evalFunctions(i,c    );
      
      % Check for function evaluation error
      if z.err == 0
        
        % Set remaining trial values
        z.evalSlacks(p,i);
        z.evalMerit (  i);
        
        % Check for nonlinear fraction-to-boundary violation
        ftb = 0;
        if i.nE > 0, ftb = ftb + sum(z.r1<min(p.ls_frac,z.mu)*r1) + sum(z.r2<min(p.ls_frac,z.mu)*r2); end;
        if i.nI > 0, ftb = ftb + sum(z.s1<min(p.ls_frac,z.mu)*s1) + sum(z.s2<min(p.ls_frac,z.mu)*s2); end;
        
        % Check Armijo condition
        if ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*max(d.qtred,0)
          
          % Reset variables, set flag, and return
          z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,z.f,z.cE,z.cI,z.phi); b = 1; return;
          
        elseif z.evalViolation(i,z.cE,z.cI) < v
          
          % Reset variables and return
          z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,f,cE,cI,phi); return;
          
        else
          
          % Reset variables (but leave constraint values for second-order correction)
          z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,f,z.cE,z.cI,phi);
          
        end
        
      else
        
        % Reset variables and return
        z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,f,cE,cI,phi); return;
        
      end
      
      % Recompute slacks for second order correction
      z.evalSlacks(p,i);
      
      % Evaluate trial primal-dual right-hand side vector
      z.evalNewtonRhs(i);
      
      % Store current direction values
      dx  = d.x ;
      dr1 = d.r1; dr2 = d.r2; dlE = d.lE;
      ds1 = d.s1; ds2 = d.s2; dlI = d.lI;
      dx_norm = d.x_norm;
      dl_norm = d.l_norm;
      
      % Evaluate search direction
      d.evalNewtonStep(i,z);
      
      % Set trial direction
      d.setDirection(i,a.p*dx+d.x,a.p*dr1+d.r1,a.p*dr2+d.r2,a.p*ds1+d.s1,a.p*ds2+d.s2,a.d*dlE+d.lE,a.d*dlI+d.lI,norm(a.p*dx+d.x),norm([a.d*dlE+d.lE;a.d*dlI+d.lI]));
      
      % Set trial point
      z.updatePoint  (i,  d,a);
      z.evalFunctions(i,c    );
      
      % Check for function evaluation error
      if z.err == 0
        
        % Set remaining trial values
        z.evalSlacks(p,i);
        z.evalMerit (  i);
        
        % Check for nonlinear fraction-to-boundary violation
        ftb = 0;
        if i.nE > 0, ftb = ftb + sum(z.r1<min(p.ls_frac,z.mu)*r1) + sum(z.r2<min(p.ls_frac,z.mu)*r2); end;
        if i.nI > 0, ftb = ftb + sum(z.s1<min(p.ls_frac,z.mu)*s1) + sum(z.s2<min(p.ls_frac,z.mu)*s2); end;
        
        % Check Armijo condition
        if ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*max(d.qtred,0)
          
          % Reset variables, set flag, and return
          z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,z.f,z.cE,z.cI,z.phi); b = 2; return;
          
        else
          
          % Reset variables
          z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,f,cE,cI,phi);
          
        end
        
      else
        
        % Reset variables
        z.setPrimals(i,x,r1,r2,s1,s2,lE,lI,f,cE,cI,phi);
        
      end
      
      % Reset direction
      d.setDirection(i,dx,dr1,dr2,ds1,ds2,dlE,dlI,dx_norm,dl_norm);
      
      % Reduce steplength
      a.p = p.ls_factor*a.p;
      
    end
    
  end
  
end
