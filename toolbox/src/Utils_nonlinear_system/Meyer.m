function RES = Meyer( what, x )
  global NF_eval NJF_eval;

  W = [ 34780, 28610, 23650, 19630, 16370, 13720, 11540, 9744, 8261, 7030, 6005, 5147, 4427, 3820, 3307, 2872 ];
  
  switch what

  case {'name','title'}
    RES = 'The Meyer Function.';

  case {'n','size'}
    RES = 3;

  case {'init','start'}
    RES = [ 0.02, 4000, 250 ]';

  case 'solution'
    RES = [];

  case 'F'
    NF_eval = NF_eval + 1;
    RES     = eval_F( x, 16, W );

  case 'JF'
    NJF_eval = NJF_eval + 1;
    RES      = eval_JF( x, 16, W );

  otherwise
    error(['The Meyer Function, bad argument:' what]);
  end;

end

%*****************************************************************************80
function f = eval_f ( x )
  global Meyer_n W
  n = Meyer_n;
  f = 0;
  for k = 1:n
    f = f + ( x(1)*exp(x(2)/(x(3)+45+5*k) )-W(k) )^2;
  end
end

%*****************************************************************************80
function g = eval_F ( x, n, W )

  g = zeros ( 3, 1 );

  for i = 1 : n
    den = x(3)+45+5*i;
    dfidx1 = exp( x(2)/den);
    
    fi =  x(1)*dfidx1 - W(i);
    
    dfidx2 = x(1)*dfidx1/den;
    dfidx3 = -x(2)*dfidx2/den;

    g(1) = g(1) + 2.0 * fi * dfidx1;
    g(2) = g(2) + 2.0 * fi * dfidx2;
    g(3) = g(3) + 2.0 * fi * dfidx3;
    
  end
end

%*****************************************************************************80
function h = eval_JF ( x, n, W )

  h = zeros ( 3, 3 );
  for i = 1 : n
    den = x(3)+45+5*i;
    esp = exp( x(2)/( x(3)+45+5*i ));

    h(1,1) = h(1,1) + 2.0 * esp^2;
    h(1,2) = h(1,2) + 2.0 * esp*(2*x(1)*esp-W(i))/den;
    h(1,3) = h(1,3) - 2.0 * x(2)*esp*(2*x(1)*esp-W(i))/den^2;

    h(2,2) = h(2,2) + 2.0 * x(1)*esp*(2*x(1)*esp-W(i))/den^2;
    h(2,3) = h(2,3) - 2.0 * x(1)*esp*(  esp*x(1)*(2*x(2)+den) - W(i)*(x(2)+den))/den^3;

    h(3,3) = h(3,3) + 2.0 * x(1)*x(2)*esp*(  2*x(1)*esp*(x(2)+den) - W(i)*(x(2)+2*den) )/den^4;
   
  end
  
  h(2,1) = h(1,2);
  h(3,1) = h(1,3);
  h(3,2) = h(3,2);
end
