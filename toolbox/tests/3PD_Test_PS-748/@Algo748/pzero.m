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
function res = pzero( ~, P, FP )
  % Uses cubic inverse interpolation of f(x) at a, b, d, and e to
  % get an approximate root of f(x).
  % This procedure is a slight modification of Aitken-Neville
  % algorithm for interpolation described by Stoer and Bulirsch
  % in "Introduction to numerical analysis" springer-verlag. new york (1980).
  %
  a = P(1); fa = FP(1);
  b = P(2); fb = FP(2);
  d = P(3); fd = FP(3);
  e = P(4); fe = FP(4);
  %Q11 = (d-e)*fd/(fe-fd);
  %Q21 = (e-d)*fb/(fd-fb);
  %Q31 = (a-b)*fa/(fb-fa);
  %D21 = (b-d)*fd/(fd-fb);
  %D31 = (a-b)*fb/(fb-fa);
  %Q22 = (D21-Q11)*fb/(fe-fb);
  %Q32 = (D31-Q21)*fa/(fb-fa);
  %D32 = (D31-Q21)*fd/(fd-fa);
  %Q33 = (D32-Q22)*fa/(fe-fa);

  Q11 = (d-e)/(fe/fd-1);
  Q21 = (e-d)/(fd/fb-1);
  Q31 = (a-b)/(fb/fa-1);
  D21 = (b-d)/(1-fb/fd);
  D31 = (a-b)/(1-fa/fb);
  Q22 = (D21-Q11)/(fe/fb-1);
  Q32 = (D31-Q21)/(fb/fa-1);
  D32 = (D31-Q21)/(1-fa/fd);
  Q33 = (D32-Q22)/(fe/fa-1);

  % Calculate the output
  res = a+(Q31+Q32+Q33);
end
