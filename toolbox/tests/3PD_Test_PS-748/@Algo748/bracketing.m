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
function [a,b,c,d,fa,fb,fc,fd] = bracketing( self, P, FP )
  % Given current enclosing interval [a,b] and a number c in (a,b):
  %
  %  a) if f(c)=0 then sets the output a=c.
  %  b) Otherwise determines the new enclosing interval:
  %     [a,b]=[a,c] or [a,b]=[c,b].
  %     also updates the termination criterion corresponding
  %     to the new enclosing interval.
  %
  % Adjust c if (b-a) is very small or if c is very close to a or b.
  a = P(1); fa = FP(1);
  b = P(2); fb = FP(2);
  c = P(3);

  tol = 0.7*self.m_tolerance;
  ba  = b - a;
  if ba <= 2*tol
    c = a+0.5*ba;
  elseif c <= a+tol
    c = a+tol;
  elseif c >= b-tol
    c = b-tol;
  end

  fc = self.f_evaluate( c );

  % If f(c)=0, then set a=c and return.
  % This will terminate the procedure.

  if fc == 0
    a = c; fa = 0;
    d = 0; fd = 0;
  else
    % If f(c) is not zero, then determine the new enclosing interval.
    if fa * fc < 0
      % D <-- B <-- C
      d  = b;  b  = c;
      fd = fb; fb = fc;
    else
      % D <-- A <-- C
      d  = a;  a  = c;
      fd = fa; fa = fc;
    end
    % update the termination criterion according to the new enclosing interval.
    if abs(fb) <= abs(fa)
      self.set_tolerance(b);
    else
      self.set_tolerance(a);
    end
  end
end
