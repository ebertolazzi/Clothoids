% Author      : Frank E. Curtis
% Description : Class for reading in AMPL model.

% Input class
classdef Input < handle
  
  % Class properties (private set access)
  properties (SetAccess = private)
    
    id % problem name
    n0 % number of original formulation variables
    I1 % indices of free variables
    I2 % indices of fixed variables
    I3 % indices of lower bounded variables
    I4 % indices of upper bounded variables
    I5 % indices of lower and upper bounded variables
    I6 % indices of equality constraints
    I7 % indices of lower bounded constraints
    I8 % indices of upper bounded constraints
    I9 % indices of lower and upper bounded constraints
    b2 % right-hand side of fixed variables
    l3 % right-hand side of lower bounded variables
    u4 % right-hand side of upper bounded variables
    l5 % right-hand side of lower half of lower and upper bounded variables
    u5 % right-hand side of upper half of lower and upper bounded variables
    b6 % right-hand side of equality constraints
    l7 % right-hand side of lower bounded constraints
    u8 % right-hand side of upper bounded constraints
    l9 % right-hand side of lower half of lower and upper bounded constraints
    u9 % right-hand side of upper half of lower and upper bounded constraints
    n1 % number of free variables
    n2 % number of fixed variables
    n3 % number of lower bounded variables
    n4 % number of upper bounded variables
    n5 % number of lower and upper bounded variables
    n6 % number of equality constraints
    n7 % number of lower bounded constraints
    n8 % number of upper bounded constraints
    n9 % number of lower and upper bounded constraints
    nV % number of variables
    nI % number of inequality constraints
    nE % number of equality constraints
    nA % size of primal-dual matrix
    x0 % initial point
    vi % counter for invalid bounds
    f_orig % original objective
    c_orig % original constraints
    g_orig % original gradient of objective
    J_orig % original jacobian of constraints
    H_orig % original hessian of the lagrangian

  end
  
  % Class methods
  methods
    
    % Constructor
    function i = Input(p, name, f_orig, c_orig, g_orig, J_orig, H_orig, x0, bl, bu, l, cl, cu)
      
      % Set problem identity
      i.id = strtrim(name);
      
      % store function pointers to original problem functions
      i.f_orig = f_orig;
      i.c_orig = c_orig;
      i.g_orig = g_orig;
      i.J_orig = J_orig;
      i.H_orig = H_orig;
      
      % Set number of original formulation variables
      i.n0 = length(x0);
      
      % Find index sets
      i.I1 = find(bl <= -p.rhs_bnd & bu >= p.rhs_bnd);
      i.I2 = find(bl ==  bu);
      i.I3 = find(bl >  -p.rhs_bnd & bu >= p.rhs_bnd);
      i.I4 = find(bl <= -p.rhs_bnd & bu <  p.rhs_bnd);
      i.I5 = find(bl >  -p.rhs_bnd & bu <  p.rhs_bnd & bl ~= bu);
      i.I6 = find(cl ==  cu);
      i.I7 = find(cl >  -p.rhs_bnd & cu >= p.rhs_bnd);
      i.I8 = find(cl <= -p.rhs_bnd & cu <  p.rhs_bnd);
      i.I9 = find(cl >  -p.rhs_bnd & cu <  p.rhs_bnd & cl ~= cu);
      
      % Set right-hand side values
      i.b2 = bl(i.I2);
      i.l3 = bl(i.I3);
      i.u4 = bu(i.I4);
      i.l5 = bl(i.I5);
      i.u5 = bu(i.I5);
      i.b6 = cl(i.I6);
      i.l7 = cl(i.I7);
      i.u8 = cu(i.I8);
      i.l9 = cl(i.I9);
      i.u9 = cu(i.I9);
      
      % Set sizes of index sets
      i.n1 = length(i.I1);
      i.n2 = length(i.I2);
      i.n3 = length(i.I3);
      i.n4 = length(i.I4);
      i.n5 = length(i.I5);
      i.n6 = length(i.I6);
      i.n7 = length(i.I7);
      i.n8 = length(i.I8);
      i.n9 = length(i.I9);
      
      % Initialize number of invalid bounds
      i.vi = 0;
      
      % Count invalid bounds
      if i.n2 > 0, i.vi = i.vi + sum(i.b2 <= -p.rhs_bnd);
                   i.vi = i.vi + sum(i.b2 >=  p.rhs_bnd); end;
      if i.n3 > 0, i.vi = i.vi + sum(i.l3 >=  p.rhs_bnd); end;
      if i.n4 > 0, i.vi = i.vi + sum(i.u4 <= -p.rhs_bnd); end;
      if i.n5 > 0, i.vi = i.vi + sum(i.l5 >=  p.rhs_bnd);
                   i.vi = i.vi + sum(i.u5 <= -p.rhs_bnd);
                   i.vi = i.vi + sum(i.l5 >   i.u5     ); end;
      if i.n6 > 0, i.vi = i.vi + sum(i.b6 <= -p.rhs_bnd);
                   i.vi = i.vi + sum(i.b6 >=  p.rhs_bnd); end;
      if i.n7 > 0, i.vi = i.vi + sum(i.l7 >=  p.rhs_bnd); end;
      if i.n8 > 0, i.vi = i.vi + sum(i.u8 <= -p.rhs_bnd); end;
      if i.n9 > 0, i.vi = i.vi + sum(i.l9 >=  p.rhs_bnd);
                   i.vi = i.vi + sum(i.u9 <= -p.rhs_bnd);
                   i.vi = i.vi + sum(i.l9 >   i.u9     ); end;
      
      % Set number of variables and constraints
      i.nV = i.n1 + i.n3 + i.n4 + i.n5;
      i.nI = i.n3 + i.n4 + 2*i.n5 + i.n7 + i.n8 + 2*i.n9;
      i.nE = i.n6;
      
      % Set size of primal-dual matrix
      i.nA = i.nV + 3*i.nE + 3*i.nI;
      
      % Set initial point
      i.x0 = [x0(i.I1); x0(i.I3); x0(i.I4); x0(i.I5);];
      
    end
    
  end
  
end
