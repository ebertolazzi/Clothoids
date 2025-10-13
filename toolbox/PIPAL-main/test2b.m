%% Rosen–Suzuki con fmincon (inline)
clear; clc; close all;

x0 = [1;1;1;1];

fun = @(x) deal( x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 ...
                 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4), ...
                 [2*x(1)-5; 2*x(2)-5; 4*x(3)-21; 2*x(4)+7] );

nonlcon = @(x) rosen_constraints(x);

options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true, ...
    'Display','iter-detailed', ...
    'MaxIterations',1000, ...
    'OptimalityTolerance',1e-12, ...
    'ConstraintTolerance',1e-12);

[x_opt,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,[],[],[],[],[],[],nonlcon,options);

fprintf('\nx* = [%g %g %g %g]\n', x_opt);
fprintf('f(x*) = %.6f\n', fval);
fprintf('Errore rispetto sol nota [0 1 2 -1]: %.3e\n', norm(x_opt-[0;1;2;-1]));


function [f, g] = rosen_suzuki_objective(x)
f = x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4);
if nargout>1
    g = [2*x(1)-5; 2*x(2)-5; 4*x(3)-21; 2*x(4)+7];
end
end


function [c, ceq, GC, GCeq] = rosen_constraints(x)
    % Disuguaglianze c <= 0
    c = [ ...
        x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(1)-x(2)+x(3)-x(4)-8;
        x(1)^2 + 2*x(2)^2 + x(3)^2 + 2*x(4)^2 - x(1) - x(4) - 10;
        2*x(1)^2 + x(2)^2 + x(3)^2 + 2*x(1) - x(2) - x(4) - 5 ];

    % Uguaglianza ceq = 0
    ceq = x(1)+x(2)+x(3)+x(4)-5;

    % Gradienti
    if nargout > 2
        GC = [ 2*x(1)+1, 2*x(1)-1, 4*x(1)+2;
               2*x(2)-1, 4*x(2),   2*x(2)-1;
               2*x(3)+1, 2*x(3),   2*x(3);
               2*x(4)-1, 4*x(4)-1, -1 ];
        GCeq = [1;1;1;1];
    end
end
