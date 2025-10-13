%% test_pipal_rosensuzuki_cle0.m
% Rosen–Suzuki riscritto con c(x) <= 0 (formato richiesto da PIPAL)
% Minimizzazione di:
% f(x) = x1^2 + x2^2 + 2 x3^2 + x4^2 - 5 x1 - 5 x2 - 21 x3 + 7 x4
% Vincoli (forma c <= 0):
% c1(x) = x1^2 + x2^2 + x3^2 + x4^2 + x1 - x2 + x3 - x4 - 8    <= 0
% c2(x) = x1^2 + 2 x2^2 + x3^2 + 2 x4^2 - x1 - x4 - 10         <= 0
% c3(x) = 2 x1^2 + x2^2 + x3^2 + 2 x1 - x2 - x4 - 5            <= 0
% c4(x) = x1 + x2 + x3 + x4 - 5                                = 0

clear; clc; close all;
addpath('PIPAL/src');   % modifica se il percorso è diverso

name = 'RosenSuzuki_cle0';

% --- obiettivo
f_orig = @(x) x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 ...
              - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4);

% gradiente obiettivo
g_orig = @(x) [2*x(1)-5; 2*x(2)-5; 4*x(3)-21; 2*x(4)+7];

% --- vincoli in forma c(x) <= 0
c_orig = @(x) [ ...
    (x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(1) - x(2) + x(3) - x(4)) - 8;
    (x(1)^2 + 2*x(2)^2 + x(3)^2 + 2*x(4)^2 - x(1) - x(4)) - 10;
    (2*x(1)^2 + x(2)^2 + x(3)^2 + 2*x(1) - x(2) - x(4)) - 5;
    (x(1) + x(2) + x(3) + x(4) - 5) ];

% Jacobiano dei vincoli (righe = grad c_i^T)
J_orig = @(x) [ ...
    2*x(1) + 1,   2*x(2) - 1,  2*x(3) + 1,  2*x(4) - 1;    % grad c1
    2*x(1) - 1,   4*x(2),      2*x(3),      4*x(4) - 1;    % grad c2
    4*x(1) + 2,   2*x(2) - 1,  2*x(3),      -1;            % grad c3
    1,            1,           1,           1 ];          % grad c4 (uguaglianza)

% Hessiano del lagrangiano: H = Hf + sum_i lambda_i * H(ci)
% Hf = diag(2,2,4,2)
% H(c1) = diag(2,2,2,2)
% H(c2) = diag(2,4,2,4)
% H(c3) = diag(4,2,2,0)
H_orig = @(x,lambda) diag([2,2,4,2]) + ...
                     lambda(1)*diag([2,2,2,2]) + ...
                     lambda(2)*diag([2,4,2,4]) + ...
                     lambda(3)*diag([4,2,2,0]);  % lambda(4) multiplies linear constraint => hess 0

%H_orig = @(x,lambda) eye(4);
H_orig = @(x,lambda) zeros(4,4);

% --- bounds, iniziale, e limiti dei vincoli (cl, cu)
x0 = [1; 1; 1; 1];               % punto iniziale (puoi cambiarlo)
bl = -inf(4,1);
bu =  inf(4,1);

% Per PIPAL: cl <= c(x) <= cu
% qui le prime 3 sono disuguaglianze c_i(x) <= 0  => cl = -inf, cu = 0
% l'ultima è uguaglianza c4(x)=0 => cl=0, cu=0
l  = [];
cl = [-inf; -inf; -inf; 0];
cu = [   0;    0;    0; 0];

outfile = 'rosensuzuki_cle0_output.txt';
algorithm = 'PIPAL_Default';

% --- istanzia e lancia
P = Pipal(name, f_orig, c_orig, g_orig, J_orig, H_orig, ...
          x0, bl, bu, l, cl, cu, outfile, algorithm );

P.optimize();

% --- risultati e diagnostica finale (senza accedere a proprietà private)
x_opt = P.getSolution();
x = x_opt.x;

fprintf('\n==== RISULTATI (c(x) <= 0) ====\n');
fprintf('x* = [%g, %g, %g, %g]\n', x(1), x(2), x(3), x(4));
fprintf('f(x*) = %.6f\n', f_orig(x));
fprintf('Errore rispetto sol nota [0 1 2 -1]: %.6e\n', norm(x - [0;1;2;-1]));

% statistiche
fval = f_orig(x);
grad = g_orig(x);
cval = c_orig(x);
J = J_orig(x);
if isfield(x_opt,'l') && ~isempty(x_opt.l)
    lambda = x_opt.l;
    kkt = norm(grad + J' * lambda);
else
    lambda = [];
    kkt = norm(grad);
end

disp('=== STATISTICHE FINALI ===');
fprintf('f(x*) = %.6f\n', fval);
fprintf('||grad|| = %.3e\n', norm(grad));
fprintf('Violazione vincoli (max(0,c)) = %.3e\n', max(0, max(cval)));
fprintf('Residuo KKT = %.3e\n', kkt);
