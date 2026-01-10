%*****************************************************************************80
%
%% Evaluates the objective function for problem 5.
%
%  Licensing:
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%    24 December 2011
%
%  Author:
%    Enrico Bertolazzi
%    Dipartimento di Ingegneria Meccanica e Strutturale
%    Universita` degli Studi di Trento
%    Via Mesiano 77, I-38123 Trento, Italy
%    email: enrico.bertolazzi@unitn.it
%
%    based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
%
%  Reference:
%    Richard Brent,
%    Algorithms for Minimization with Derivatives,
%    Dover, 2002,
%    ISBN: 0-486-41998-3,
%    LC: QA402.5.B74.
%
function RES = Minimal11(what,x)
  global NF_eval NJF_eval;

  switch what

  case {'name','title'}
    RES = 'The Minimal function n7.';

  case {'n','size'}
    RES = 11;

  case {'init','start'}
    RES = ones(11,1);

  case 'solution'
    RES = 1.000009197 * ones(11,1);

  case 'F'
    NF_eval = NF_eval + 1;
    RES     = eval_F( x, 11 );

  case 'JF'
    NJF_eval = NJF_eval + 1;
    RES      = eval_JF( x, 11 );

  otherwise
    error(['The Minimal function n7:' what]);
  end;

end

%*****************************************************************************80
function f = eval_f ( x, n )
  f = 0.0;
  for i = 1 : n
    ll = log(x(i));
    ee = exp(x(i));
    fi =  ll + ee - sqrt((ll - ee)^2 + 1e-4);
    f = f + fi * fi;
  end
end

%*****************************************************************************80
function g = eval_F ( x, n )
  g = zeros ( n, 1 );
  for i = 1 : n
    ll = log(x(i));
    ee = exp(x(i));
    dll = 1/x(i);
    fi =  ll + ee - sqrt((ll - ee)^2 + 1e-4);
    g(i) = 2 * fi * (dll + ee - 1/sqrt((ll + ee)^2 + 1e-4)*((ll - ee)*(dll - ee)));
  end
end

%*****************************************************************************80
function h = eval_JF ( x, n )
  h = zeros( n, n );
  for i = 1 : n
    ll = log(x(i));
    ee = exp(x(i));
    h(i,i) = 2 * (1 / x(i) + ee - ((ll - ee) ^ 2 + 1e-4) ^ (-1 / 2) * (ll - ee) * (1 / x(i) - ee)) ^ 2 + 2 * (ll + ee - sqrt((ll - ee) ^ 2 + 1e-4)) * (-1 / x(i) ^ 2 + ee + ((ll - ee) ^ 2 + 1e-4) ^ (-3 / 2) * (ll - ee) ^ 2 * (1 / x(i) - ee) ^ 2 - ((ll - ee) ^ 2 + 1e-4) ^ (-1 / 2) * (1 / x(i) - ee) ^ 2 - ((ll - ee) ^ 2 + 1e-4) ^ (-1 / 2) * (ll - ee) * (-1 / x(i) ^ 2 - ee));
  end
end
