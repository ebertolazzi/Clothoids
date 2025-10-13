% Author      : Frank E. Curtis
% Description : Class for algorithm counters.

% Counters class
classdef Counter < handle
  
  % Class properties (private set access)
  properties (SetAccess = private)
    
    f = 0; % Function evaluation counter
    g = 0; % Gradient evaluation counter
    H = 0; % Hessian evaluation counter
    k = 0; % Iteration counter
    M = 0; % Matrix factorization counter
    
  end
  
  % Class methods
  methods
    
    % Matrix factorization counter incrementor
    function incrementFactorizationCount(c)
      c.M = c.M + 1;
    end
    
    % Function evaluation counter incrementor
    function incrementFunctionCount(c)
      c.f = c.f + 1;
    end
    
    % Gradient evaluation counter incrementor
    function incrementGradientCount(c)
      c.g = c.g + 1;
    end
    
    % Hessian evaluation counter incrementor
    function incrementHessianCount(c)
      c.H = c.H + 1;
    end
    
    % Iteration counter incrementor
    function incrementIterationCount(c)
      c.k = c.k + 1;
    end
    
  end
  
end