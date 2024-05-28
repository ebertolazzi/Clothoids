%-------------------------------------------------------+
%                                                       |
% Copyright (C) 2022                                    |
%                                                       |
%        , __                 , __                      |
%       /|/  \               /|/  \                     |
%        | __/ _   ,_         | __/ _   ,_              |
%        |   \|/  /  |  |   | |   \|/  /  |  |   |      |
%        |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/     |
%                          /|                   /|      |
%                          \|                   \|      |
%                                                       |
%     Enrico Bertolazzi                                 |
%     Dipartimento di Ingegneria Industriale            |
%     Universita` degli Studi di Trento                 |
%     email: enrico.bertolazzi@unitn.it                 |
%                                                       |
%-------------------------------------------------------+
function [sol,iter] = eval( self, a, b, fun )
  self.m_num_fun_eval  = 0;
  self.m_num_iter_done = 0;
  self.m_fun           = fun;
  % check for trivial solution
  fa = self.f_evaluate( a ); if fa == 0; sol = a; iter = 0; return; end
  fb = self.f_evaluate( b ); if fb == 0; sol = b; iter = 0; return; end
  % check if solution can exists
  if fa*fb > 0
    iter = -1;
    sol  = 0;
    return;
  end
  % Finds either an exact solution or an approximate solution
  % of the equation f(x)=0 in the interval [a,b].
  %
  % At the beginning of each iteration, the current enclosing interval
  % is recorded as [a0,b0].
  % The first iteration is simply a secant step.
  % Starting with the second iteration, three steps are taken in each iteration.
  % First two steps are either quadratic interpolation
  % or cubic inverse interpolation.
  % The third step is a double-size secant step.
  % If the diameter of the enclosing interval obtained after
  % those three steps is larger than 0.5*(b0-a0),
  % then an additional bisection step will be taken.
  e = NaN; fe = NaN; % Dumb values
  % Until f(left) or f(right) are infinite perform bisection
  while ~( isfinite(fa) && isfinite(fb) )
    self.m_num_iter_done = self.m_num_iter_done+1;
    c  = (a+b)/2;
    fc = self.f_evaluate( c );
    if fc == 0
      sol  = c;
      iter = self.m_num_iter_done;
      return;
    end
    if fa*fc < 0
      % --> [a,c]
      b = c; fb = fc;
    else
      % --> [c,b]
      a = c; fa = fc;
    end
    if ( abs(fb) <= abs(fa) )
      sol = b;
      self.set_tolerance(b);
    else
      sol = a;
      self.set_tolerance(a);
    end
    if (b-a) <= self.m_tolerance
      iter = self.m_num_iter_done;
      return; % found solution
    end
  end
  ba  = b-a;
  fba = fb-fa;
  R   = ba/fba;
  if abs(fb) < abs(fa)
    c = b+fb*R;
  else
    c = a-fa*R;
  end
  % impedisce m_c troppo vicino ad m_a o m_b
  c = max(min(c,b-0.1*ba),a+0.1*ba);
  %
  % Call "bracketing" to get a shrinked enclosing interval as
  % well as to update the termination criterion.
  % Stop the procedure if the criterion is satisfied or the
  % exact solution is obtained.
  %
  [a,b,c,d,fa,fb,fc,fd] = self.bracketing( [a,b,c], [fa,fb] );
  if fa == 0
    iter = self.m_num_iter_done;
    sol  = a;
    return;
  end
  % Iteration starts.
  % The enclosing interval before executing the iteration is recorded as [a0, b0].
  converged = false;
  while ~converged
    self.m_num_iter_done = self.m_num_iter_done+1;
    A0   = a;
    B0   = b;
    % Calculates the termination criterion.
    % Stops the procedure if the criterion is satisfied.
    if abs(fb) <= abs(fa)
      self.set_tolerance(b);
    else
      self.set_tolerance(a);
    end
    converged = (b-a) <= self.m_tolerance;
    if converged; break; end
    %
    % Starting with the second iteration, in the first two steps, either
    % quadratic interpolation is used by calling the subroutine "newtonquadratic"
    % or the cubic inverse interpolation is used by calling the subroutine
    % "pzero". in the following, if "prof" is not equal to 0, then the
    % four function values "fa", "fb", "fd", and "fe" are distinct, and
    % hence "pzero" will be called.
    %
    do_newton_quadratic = false;
    if isnan(fe)
      do_newton_quadratic = true;
    elseif ~self.all_different([fa,fb,fd,fe])
      do_newton_quadratic = true;
    else
      c = self.pzero([a,b,d,e],[fa,fb,fd,fe]);
      if (c-a)*(c-b) >= 0
        do_newton_quadratic = true;
      end
    end
    if do_newton_quadratic
      c = self.newton_quadratic(2,[a,b,d],[fa,fb,fd]);
    end
    e  = d;
    fe = fd;
    %
    % Call subroutine "bracketing" to get a shrinked enclosing interval as
    % well as to update the termination criterion. stop the procedure
    % if the criterion is satisfied or the exact solution is obtained.
    %
    [a,b,c,d,fa,fb,fc,fd] = self.bracketing( [a,b,c], [fa,fb] );
    converged = fa == 0 || (b-a) <= self.m_tolerance;
    if converged; break; end
    do_newton_quadratic = false;
    if ~self.all_different([fa,fb,fd,fe])
      do_newton_quadratic = true;
    else
      c = self.pzero([a,b,d,e],[fa,fb,fd,fe]);
      if (c-a)*(c-b) >= 0
        do_newton_quadratic = true;
      end
    end
    if do_newton_quadratic
      c = self.newton_quadratic(3,[a,b,d],[fa,fb,fd]);
    end
    %
    % Call subroutine "bracketing" to get a shrinked enclosing interval as
    % well as to update the termination criterion. stop the procedure
    % if the criterion is satisfied or the exact solution is obtained.
    %
    [a,b,c,d,fa,fb,fc,fd] = self.bracketing( [a,b,c], [fa,fb] );
    converged = fa == 0 || (b-a) <= self.m_tolerance;
    if converged; break; end
    e  = d;
    fe = fd;
    % Takes the double-size secant step.
    if abs(fa) < abs(fb)
      u = a; fu = fa;
    else
      u = b; fu = fb;
    end
    c = u-2*(fu/(fb-fa))*(b-a);
    if abs(c-u) > 0.5*(b-a)
      c = a+0.5*(b-a);
    end
    %
    % Call subroutine "bracketing" to get a shrinked enclosing interval as
    % well as to update the termination criterion. stop the procedure
    % if the criterion is satisfied or the exact solution is obtained.
    %
    [a,b,c,d,fa,fb,fc,fd] = self.bracketing( [a,b,c], [fa,fb] );
    converged = fa == 0 || (b-a) <= self.m_tolerance;
    if converged; break; end
    %
    % Determines whether an additional bisection step is needed.
    % Takes it if necessary.
    %
    if (b-a) < self.m_mu*(B0-A0)
      continue;
    end
    e  = d;
    fe = fd;
    %
    % Call subroutine "bracketing" to get a shrinked enclosing interval as
    % well as to update the termination criterion. stop the procedure
    % if the criterion is satisfied or the exact solution is obtained.
    %
    c = a+0.5*(b-a);
    [a,b,c,d,fa,fb,fc,fd] = self.bracketing( [a,b,c], [fa,fb] );
    converged = fa == 0 || (b-a) <= self.m_tolerance;
  end
  % terminates the procedure and return the "root".
  sol = a;
  if converged
    iter = self.m_num_iter_done;
  else
    iter = -self.m_num_iter_done;
  end
end
