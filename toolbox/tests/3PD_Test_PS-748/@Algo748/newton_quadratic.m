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
function res = newton_quadratic( ~, niter, P, FP )
  % Uses `niter` newton steps to approximate the zero in (a,b) of the
  % quadratic polynomial interpolating f(x) at a, b, and d.
  % Safeguard is used to avoid overflow.
  a = P(1); fa = FP(1);
  b = P(2); fb = FP(2);
  d = P(3); fd = FP(3);
  A0 = fa;
  A1 = (fb-fa)/(b-a);
  A2 = ((fd-fb)/(d-b)-A1)/(d-a);
  % Safeguard to avoid overflow.
  if A2 == 0
    res = a-A0/A1;
  else
    % Determine the starting point of newton steps.
    if A2*fa > 0
      c = a;
    else
      c = b;
    end
    % Start the safeguarded newton steps.
    ok = true;
    for i=1:niter
      PC  = A0+(A1+A2*(c-b))*(c-a);
      PDC = A1+A2*((2*c)-(a+b));
      ok = PDC ~= 0;
      if ~ok; break; end
      c = c - PC/PDC;
    end
    if ok
      res = c;
    else
      res = a-A0/A1;
    end
  end
end
