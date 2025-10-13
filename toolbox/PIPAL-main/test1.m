%% test_rosenbrock_box.m
% Test del solver PIPAL sul problema di Rosenbrock
% vincolato al quadrato [-1,1] x [-1,1]

clear; clc; close all;

addpath('PIPAL/src');

%% test_pipal_rosenbrock.m
% Esempio d'uso del solver PIPAL (Frank E. Curtis)
% Minimizzazione della funzione di Rosenbrock nel quadrato [-1,1]^2
%% === 2. Definizione del problema ===

name = 'RosenbrockBox';

% Funzione obiettivo
f_orig = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;

% Gradiente
g_orig = @(x) [ -400*x(1)*(x(2)-x(1)^2) - 2*(1 - x(1));
                 200*(x(2) - x(1)^2) ];

% Vincoli non lineari (nessuno)
c_orig = @(x) [];

% Jacobiano dei vincoli
J_orig = @(x) [];

% Hessiano del lagrangiano (solo parte obiettivo)
H_orig = @(x,lambda) [1200*x(1)^2 - 400*x(2) + 2,   -400*x(1);
                      -400*x(1),                     200];

% Punto iniziale e limiti
x0 = [0; 0];
bl = [-1; -1];
bu = [ 1;  1];

% Nessun vincolo lineare
l  = [];
cl = [];
cu = [];

% File di output del log
outfile = 'rosenbrock_output.txt';

% Algoritmo interno
algorithm = 'PIPAL_Default';

%% === 3. Istanzia l'oggetto solver ===
P = Pipal(name, f_orig, c_orig, g_orig, J_orig, H_orig, ...
          x0, bl, bu, l, cl, cu, outfile, algorithm);

%% === 4. Esegui l'ottimizzazione ===
P.optimize();

%% === 5. Recupera e mostra i risultati ===
x_opt = P.getSolution();
fprintf('\n==== RISULTATI ====\n');
fprintf('x* = [%g, %g]\n', x_opt.x(1), x_opt.x(2));
fprintf('f(x*) = %.6f\n', f_orig(x_opt.x));
