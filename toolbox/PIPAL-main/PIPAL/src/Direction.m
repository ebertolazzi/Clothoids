% Author      : Frank E. Curtis
% Description : Class for direction quantities and routines; methods for
%               evaluating search directions and model values/reductions.

% Direction class
classdef Direction < handle
  
  % Class properties (private set access)
  properties (SetAccess = private)
    
    x       % Primal direction
    x_norm  % Primal direction norm value
    x_norm_ % Primal direction norm last value
    r1      % Equality constraint slack direction
    r2      % Equality constraint slack direction
    lE      % Equality constraint multiplier direction
    s1      % Inequality constraint slack direction
    s2      % Inequality constraint slack direction
    lI      % Inequality constraint multiplier direction
    l_norm  % Constraint multiplier direction norm
    lred0   % Penalty-interior-point linear model value for zero penalty parameter
    ltred0  % Penalty-interior-point linear model reduction value for zero penalty parameter
    ltred   % Penalty-interior-point linear model reduction value
    qtred   % Penalty-interior-point quadratic model reduction value
    m       % Quality function value
    
  end
  
  % Class methods
  methods
    
    % Constructor
    function d = Direction
      
      % Initialize last direction norm
      d.x_norm_ = inf;
      
    end
    
    % Evaluate linear combination of directions
    function evalLinearCombination(d,i,d1,d2,d3,a)
      
      % Evaluate linear combinations
      d.x  = a(1)*d1.x  + a(2)*d2.x  + a(3)*d3.x;
      if i.nE > 0, d.r1 = a(1)*d1.r1 + a(2)*d2.r1 + a(3)*d3.r1;
                   d.r2 = a(1)*d1.r2 + a(2)*d2.r2 + a(3)*d3.r2; end;
      if i.nI > 0, d.s1 = a(1)*d1.s1 + a(2)*d2.s1 + a(3)*d3.s1;
                   d.s2 = a(1)*d1.s2 + a(2)*d2.s2 + a(3)*d3.s2; end;
      if i.nE > 0, d.lE = a(1)*d1.lE + a(2)*d2.lE + a(3)*d3.lE; end;
      if i.nI > 0, d.lI = a(1)*d1.lI + a(2)*d2.lI + a(3)*d3.lI; end;
      
      % Evaluate primal direction norm
      d.x_norm = norm(d.x);
      
      % Evaluate dual direction norm
      d.l_norm = norm([d.lE;d.lI]);
      
    end
    
    % Evaluate model and model reductions
    function evalModels(d,i,z)
      
      % Evaluate reduction in linear model of penalty-interior-point objective for zero penalty parameter
      d.lred0 = 0;
      if i.nE > 0, d.lred0 = d.lred0 - sum([1-z.mu./z.r1; 1-z.mu./z.r2].*[d.r1; d.r2]); end;
      if i.nI > 0, d.lred0 = d.lred0 - sum([0-z.mu./z.s1; 1-z.mu./z.s2].*[d.s1; d.s2]); end;
      
      % Evaluate remaining quantities only for nonzero penalty parameter
      if z.rho > 0
        
        % Evaluate reduction in linear model of merit function for zero penalty parameter
        d.ltred0 = 0;
        if i.nE > 0, d.ltred0 = d.ltred0 - (1/2)*full(sum(((1-z.mu./z.r1).*(-1+z.cE./(sqrt(z.cE.^2 +   z.mu^2)))+(1-z.mu./z.r2).*(1+z.cE./(sqrt(z.cE.^2 +   z.mu^2)))).*(z.JE*d.x))); end;
        if i.nI > 0, d.ltred0 = d.ltred0 - (1/2)*full(sum(((0-z.mu./z.s1).*(-1+z.cI./(sqrt(z.cI.^2 + 4*z.mu^2)))+(1-z.mu./z.s2).*(1+z.cI./(sqrt(z.cI.^2 + 4*z.mu^2)))).*(z.JI*d.x))); end;
        
        % Evaluate reduction in linear model of merit function
        d.ltred = -z.rho*z.g'*d.x + d.ltred0;
        
        % Evaluate reduction in quadratic model of merit function
        d.qtred = d.ltred - (1/2)*d.x'*z.H*d.x;
        if i.nE > 0, Jd = z.JE*d.x; Dinv = z.r1./(1+z.lE)+z.r2./(1-z.lE); d.qtred = d.qtred - (1/2)*Jd'*(Jd./Dinv); end;
        if i.nI > 0, Jd = z.JI*d.x; Dinv = z.s1./(0+z.lI)+z.s2./(1-z.lI); d.qtred = d.qtred - (1/2)*Jd'*(Jd./Dinv); end;
        
        % Initialize quality function vector
        vec = zeros(i.nV+2*i.nE+2*i.nI,1);
        
        % Set gradient of objective
        vec(1:i.nV) = z.rho*z.g;
        
        % Set gradient of Lagrangian for constraints
        if i.nE > 0, vec(1:i.nV) = vec(1:i.nV) + ((z.lE+d.lE)'*z.JE)'; end;
        if i.nI > 0, vec(1:i.nV) = vec(1:i.nV) + ((z.lI+d.lI)'*z.JI)'; end;
        
        % Set complementarity for constraint slacks
        if i.nE > 0, vec(1+i.nV       :i.nV+2*i.nE       ) = [(z.r1+d.r1).*(1 + (z.lE+d.lE)); (z.r2+d.r2).*(1 - (z.lE+d.lE))]; end;
        if i.nI > 0, vec(1+i.nV+2*i.nE:i.nV+2*i.nE+2*i.nI) = [(z.s1+d.s1).*(0 + (z.lI+d.lI)); (z.s2+d.s2).*(1 - (z.lI+d.lI))]; end;
        
        % Evaluate quality function
        d.m = norm(vec,inf);
        
      end
      
    end
    
    % Evaluate Newton step
    function evalNewtonStep(d,i,z)
      
      % Evaluate direction
      dir = z.AS(:,z.AP)*(z.AL'\(z.AD\(z.AL\(z.AS(z.AP,:)*(-z.b)))));
      
      % Parse direction
      d.x               = dir(1                              :i.nV                              );
      if i.nE > 0, d.r1 = dir(1+i.nV                         :i.nV+i.nE                         );
                   d.r2 = dir(1+i.nV+i.nE                    :i.nV+i.nE+i.nE                    ); end;
      if i.nI > 0, d.s1 = dir(1+i.nV+i.nE+i.nE               :i.nV+i.nE+i.nE+i.nI               );
                   d.s2 = dir(1+i.nV+i.nE+i.nE+i.nI          :i.nV+i.nE+i.nE+i.nI+i.nI          ); end;
      if i.nE > 0, d.lE = dir(1+i.nV+i.nE+i.nE+i.nI+i.nI     :i.nV+i.nE+i.nE+i.nI+i.nI+i.nE     ); end;
      if i.nI > 0, d.lI = dir(1+i.nV+i.nE+i.nE+i.nI+i.nI+i.nE:i.nV+i.nE+i.nE+i.nI+i.nI+i.nE+i.nI); end;
      
      % Evaluate primal direction norm
      d.x_norm = norm(d.x);
      
      % Evaluate dual direction norm
      d.l_norm = norm([d.lE;d.lI]);
      
    end
    
    % Evaluate search direction quantities
    function evalStep(d,p,i,c,z,a)
      
      % Reset maximum exponent for interior-point parameter increases
      p.resetMuMaxExp;
      
      % Update penalty-interior-point parameters based on KKT errors
      z.updateParameters(p,i);
      
      % Evaluate matrices
      z.evalMatrices(p,i,c);
      
      % Set last penalty parameter
      z.setRhoLast(z.rho);
      
      % Check for aggressive algorithm
      if p.algorithm == 1
        
        % Check KKT memory for potential mu increase limit
        if z.kkt(2) > max(z.kkt_), p.setMuMaxExpZero; end;
        
        % Store current penalty and interior-point parameters
        rho_curr = z.rho; mu_curr = z.mu;
        
        % Evaluate trial steps
        [d1,d2,d3] = d.evalTrialSteps(i,z);
        
        % Set trial interior-point parameter values
        Mu = max(p.mu_min,min(p.mu_factor.^([p.mu_trials-1:-1:0]-p.mu_max_exp)*mu_curr,p.mu_max));
        
        % Initialize feasibility direction data
        lred0_0_mu = zeros(1,p.mu_trials);
        
        % Loop through interior-point parameter values
        for j = 1:p.mu_trials
          
          % Set penalty and interior-point parameters
          z.setRho(0); z.setMu(Mu(j));
          
          % Evaluate direction
          d.evalLinearCombination(i,d1,d2,d3,[(z.rho/rho_curr+z.mu/mu_curr-1),(1-z.mu/mu_curr),(1-z.rho/rho_curr)]);
          
          % Cut length
          d.x = min(d.x_norm_/max(d.x_norm,1),1)*d.x;
          
          % Run fraction-to-boundary
          a.fractionToBoundary(p,i,z,d);
          
          % Cut length
          d.evalTrialStepCut(i,a);
          
          % Evaluate models
          d.evalModels(i,z);
          
          % Set feasibility direction data
          lred0_0_mu(j) = d.lred0;
          
        end
        
        % Initialize updating data
        ltred0_rho_mu = zeros(p.mu_trials,1);
        qtred_rho_mu  = zeros(p.mu_trials,1);
        m_rho_mu      = zeros(p.mu_trials,1);
        
        % Initialize check
        check = 0;
        
        % Loop through penalty parameter values
        for k = 1:p.rho_trials
          
          % Set penalty parameter
          z.setRho(max(p.rho_min,(p.rho_factor^(k-1))*rho_curr));
          
          % Set last penalty parameter
          if rho_curr > z.kkt(1)^2, z.setRhoLast(z.rho); end;
          
          % Loop through interior-point parameter values
          for j = 1:p.mu_trials
            
            % Set interior-point parameter
            z.setMu(Mu(j));
            
            % Evaluate direction
            d.evalLinearCombination(i,d1,d2,d3,[(z.rho/rho_curr+z.mu/mu_curr-1),(1-z.mu/mu_curr),(1-z.rho/rho_curr)]);
            
            % Run fraction-to-boundary
            a.fractionToBoundary(p,i,z,d);
            
            % Cut steps
            d.evalTrialStepCut(i,a);
            
            % Evaluate models
            d.evalModels(i,z);
            
            % Set updating data
            ltred0_rho_mu(j) = d.ltred0;
            qtred_rho_mu(j)  = d.qtred;
            m_rho_mu(j)      = d.m;
            
            % Check updating conditions for infeasible points
            if z.v >  p.opt_err_tol && (ltred0_rho_mu(j) < p.update_con_1*lred0_0_mu(j) || qtred_rho_mu(j) < p.update_con_2*lred0_0_mu(j) || z.rho > z.kkt(1)^2), m_rho_mu(j) = inf; end;
            
            % Check updating conditions for feasible points
            if z.v <= p.opt_err_tol && qtred_rho_mu(j) < 0, m_rho_mu(j) = inf; end;
            
          end
          
          % Find minimum m for current rho
          m_min = min(m_rho_mu);
          
          % Check for finite minimum
          if m_min < inf
            
            % Loop through mu values
            for j = 1:p.mu_trials
              
              % Set mu
              mu = Mu(j);
              
              % Check condition
              if m_rho_mu(j) <= p.update_con_3*m_min, z.setMu(mu); end;
              
            end
            
            % Set condition check
            check = 1;
            
            % Break loop
            break;
            
          end
          
        end
        
        % Check conditions
        if check == 0, z.setRho(rho_curr); z.setMu(mu_curr); end;
        
        % Evaluate merit
        z.evalMerit(i);
        
      end
      
      % Evaluate primal-dual right-hand side vector
      z.evalNewtonRhs(i);
      
      % Evaluate search direction
      d.evalNewtonStep(i,z);
      
      % Evaluate models
      d.evalModels(i,z);
      
      % Store last direction norm
      d.x_norm_ = d.x_norm;
      
    end
    
    % Evaluate and store trial step
    function v = evalTrialStep(d,i)
      
      % Set direction components
                   v.x  = d.x ;
      if i.nE > 0, v.r1 = d.r1;
                   v.r2 = d.r2; end;
      if i.nI > 0, v.s1 = d.s1;
                   v.s2 = d.s2; end;
      if i.nE > 0, v.lE = d.lE; end;
      if i.nI > 0, v.lI = d.lI; end;
      
    end
    
    % Evaluate trial step cut by fraction-to-boundary rule
    function evalTrialStepCut(d,i,a)
      
      % Set direction components
                   d.x  = a.p*d.x ;
      if i.nE > 0, d.r1 = a.p*d.r1;
                   d.r2 = a.p*d.r2; end;
      if i.nI > 0, d.s1 = a.p*d.s1;
                   d.s2 = a.p*d.s2; end;
      if i.nE > 0, d.lE = a.d*d.lE; end;
      if i.nI > 0, d.lI = a.d*d.lI; end;
      
    end
    
    % Evaluate and store directions for parameter combinations
    function [d1,d2,d3] = evalTrialSteps(d,i,z)
      
      % Store current penalty and interior-point parameters
      rho_curr = z.rho;
      mu_curr  = z.mu;
      
      % Evaluate direction for current penalty and interior-point parameters
      z.setRho(rho_curr);
      z.setMu(mu_curr);
      z.evalNewtonRhs(i);
      d.evalNewtonStep(i,z);
      d1 = d.evalTrialStep(i);
      
      % Evaluate direction for zero interior-point parameter
      z.setRho(rho_curr);
      z.setMu(0);
      z.evalNewtonRhs(i);
      d.evalNewtonStep(i,z);
      d2 = d.evalTrialStep(i);
      
      % Evaluate direction for zero penalty parameter
      z.setRho(0);
      z.setMu(mu_curr);
      z.evalNewtonRhs(i);
      d.evalNewtonStep(i,z);
      d3 = d.evalTrialStep(i);
      
    end
    
    % Set direction
    function setDirection(d,i,dx,dr1,dr2,ds1,ds2,dlE,dlI,dx_norm,dl_norm)
      
      % Set primal variables
      d.x = dx;
      if i.nE > 0, d.r1 = dr1; d.r2 = dr2; d.lE = dlE; end;
      if i.nI > 0, d.s1 = ds1; d.s2 = ds2; d.lI = dlI; end;
      d.x_norm = dx_norm;
      d.l_norm = dl_norm;
      
    end
    
  end
  
end
