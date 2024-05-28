%=========================================================================%
%                                                                         %
%  Autors: Enrico Bertolazzi                                              %
%          Department of Industrial Engineering                           %
%          University of Trento                                           %
%          enrico.bertolazzi@unitn.it                                     %
%          m.fregox@gmail.com                                             %
%                                                                         %
%=========================================================================%
% Driver test program to check Clothoids lib                              %
%=========================================================================%

close all;
clc;
clear all;
format long;

% Set LaTeX as default interpreter for axis labels, ticks and legends
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',18)
set(0,'DefaultLegendFontSize',18)

%addpath('../lib');
%addpath('../bin');
%addpath('../tests');
% addpath('./Algo748')

figure_flag = false;
save_flag   = false;


DB_A = Dubins();
DB_B = Dubins();

k_max  = 0.9;
d      = 3;
x0     = -d;
y0     = 0;
xM     = 0;
yM     = 0;
xf     = d;
yf     = 0;
theta0 = -pi/2.0*0-pi;
thetaf = -pi/2.0*0-pi;
len    = [];
DlenA  = [];
DlenB  = [];
kind   = [];
epsilon = 1e-4;

max_fig = 100;
if figure_flag
  figure();
end

% createa vector of combination of
theta0_v    = linspace(-pi,pi,100);
thetaf_v    = linspace(-pi,pi,100);
k_max_v     = linspace(0.1,2,10);
theta_int_v = linspace(0,2*pi,10);

mesh        = combvec(theta0_v,thetaf_v,k_max_v,theta_int_v);

% results     = [mesh;zeros(4,size(mesh,2))];

results2 = zeros(7,size(mesh,2));

for id=1:length(mesh)
  try
    results2(:,id) = sim(mesh,id,figure_flag,save_flag);
  catch
    fprintf('Error in test %d\n',id);
    results2(:,id) = results2(:,id)*0 -1;
    fprintf('%g %g %g %g\n',mesh(:,id));
    % pause(10)
  end


end

results     = [mesh;results2];

res_table = array2table(results',...
    'VariableNames',{'theta0','thetaf','k_max','gamma','thetaM','n_iter','L','DL','typeA','typeB','time'});

% print a bit of statistics of the results
fprintf('------------------------------------------------------------\n');
fprintf('Statistics\n');
fprintf('mean n_iter = %g\n',mean(results(6,:)));
fprintf('max  n_iter = %g\n',max(results(6,:)));
fprintf('min  n_iter = %g\n',min(results(6,:)));
fprintf('STD  n_iter = %g\n',std(results(6,:)));
fprintf('mean time   = %g\n',mean(results(11,:)));
fprintf('max  time   = %g\n',max(results(11,:)));
fprintf('min  time   = %g\n',min(results(11,:)));
fprintf('STD  time   = %g\n',std(results(11,:)));
fprintf('------------------------------------------------------------\n');


% save a csv file named with date and hour with the results and a header for the names of the columns
file_name = sprintf('./img/test_Dubins75_%s.csv',datestr(now,'yyyy_mm_dd_HH_MM_SS'));

header = {'theta0','thetaf','k_max','gamma','thetaM','n_iter','L','DL','typeA','typeB','time'};
fid = fopen(file_name,'w');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',header{:});
fprintf(fid,'%g,%g,%g,%g,%g,%d,%g,%g,%d,%d,%d\n',results);
fclose(fid);

save('img/DubinsPSC.mat', 'res_table')

%%
wc = find(results(6,:)>10000);
for id = wc
  sim(mesh,id,true,false)
  pause(0.001)
end



%%
function len = len_Dubins( x0, y0, theta0, xf, yf, thetaf, k_max )
  DB = Dubins();
  DB.build( x0, y0, theta0, xf, yf, thetaf, k_max );
  [len,~,~] = DB.length();
end

%%
function DL = len_Dubins_DL( x0, y0, theta0, xf, yf, thetaf, k_max )
  DB = Dubins();
  DB.build( x0, y0, theta0, xf, yf, thetaf, k_max );
  [~,DL,~] = DB.length();
end

%%
function DR = len_Dubins_DR( x0, y0, theta0, xf, yf, thetaf, k_max )
  DB = Dubins();
  DB.build( x0, y0, theta0, xf, yf, thetaf, k_max );
  [~,~,DR] = DB.length();
end

%%
function type = type_Dubins( x0, y0, theta0, xf, yf, thetaf, k_max )
  DB = Dubins();
  DB.build( x0, y0, theta0, xf, yf, thetaf, k_max );
  type = DB.curve_type();
end

%%


% function to perform a pattern search to find the optimal thetaM0
function [root, iter] = pattern_search( func, tol, max_iter )
  % Initialize variables
  theta_max       = pi;
  theta_min       = -pi;
  theta_candidate = 0;
  numpts          = 16;
  delta_theta     = (theta_max - theta_min) / numpts;
  iter = 0;
  while ( delta_theta > tol && iter < max_iter)
    thetas = linspace(theta_min, theta_max, numpts);
    residuals = zeros(1, numpts);
    for i = 1:length(thetas)
      residuals(i) = func(thetas(i));
    end
    [~, min_idx] = min(residuals);
    theta_candidate = thetas(min_idx);
    theta_min       = theta_candidate-delta_theta;
    theta_max       = theta_candidate+delta_theta;
    delta_theta     = (theta_max - theta_min) / numpts;
    iter = iter + numpts;
  end
  root = theta_candidate;
  fprintf("Pattern: Solution found at iter %d\n",iter);
end


function [root, iter] = bisection_method(func, a, b, tol, max_iter)
  if sign(func(a)) == sign(func(b))
    error('Function has same sign at endpoints. Bisection method cannot be applied.');
  end
  iter = 0;
  while (b - a) / 2 > tol && iter < max_iter
      c = (a + b) / 2;
      if func(c) == 0
        root = c;
        return;
      end
      if sign(func(c)) == sign(func(a))
        a = c;
      else
        b = c;
      end
      iter = iter + 1;
  end
  root = (a + b) / 2;
  fprintf("Bisection: Solution found at iter %d\n",iter);
end


function xn = EXPEXP(a, F, eps, n)
  a1 = a;
  g = (F(a1 + F(a1)) - F(a1)) / F(a1);
  y = a1 * exp(-F(a1) / (a1 * g));
  h = (F(y) - F(a1)) / (y - a1);
  xn = y * exp(-F(y) / (y * h));
  i = 1;
  while F(xn) ~= 0 && F(y) ~= 0 && xn - a1 >= eps && i < n
      a1 = xn;
      i = i + 1;
      y = a1 * exp(-F(a1) / (a1 * g));
      h = (F(y) - F(a1)) / (y - al);
      xn = y * exp(-F(y) / (y * h));
      fprintf('Iteration %d: x = %g\n', i, xn);
  end
end

% Mixed strategy. Pattern search untill we have same type of solution at the border of the interval of the pattern search and the gradient changes sign. Then we use bisection method to find the root.
%

function [root, tot_iter] = mixed_strategy(func, d_func, ftype, tol, max_iter)
  % Initialize variables
  theta_max       = pi;
  theta_min       = -pi;
  theta_candidate = 0;
  numpts          = 16;
  delta_theta     = (theta_max - theta_min) / numpts;
  iter = 0;
  subiter = 0;
  while (delta_theta > tol && iter < max_iter)
    thetas = linspace(theta_min, theta_max, numpts);
    residuals = zeros(1, numpts);
    residuals = arrayfun(func, thetas);
    [~, min_idx] = min(residuals);
    theta_candidate = thetas(min_idx);
    theta_min       = theta_candidate-delta_theta;
    theta_max       = theta_candidate+delta_theta;
    delta_theta     = (theta_max - theta_min) / numpts;
    iter = iter + numpts;
    numpts = max(4,numpts/2);
    type_min   = ftype(theta_min);
    type_max   = ftype(theta_max);
    d_func_min = d_func(theta_min);
    d_func_max = d_func(theta_max);
    if (type_min(1) == type_max(1) && type_min(2) == type_max(2) && sign(d_func_min) ~= sign(d_func_max))
      [root, subiter] = bisection_method(d_func, theta_min, theta_max, tol, max_iter);
      tot_iter = iter + subiter;
      %fprintf("Mixed Strategy (pattern + bisec): Solution found at iter %d\n",tot_iter);
      return;
    end
  end
  root = theta_candidate;
  tot_iter = iter + subiter;
  %fprintf("Mixed Strategy (pattern): Solution found at iter %d\n",tot_iter);
end

function [root, tot_iter] = mixed_strategy748(func, d_func, ftype, tol, max_iter)
  % Initialize variables
  theta_max       = pi;
  theta_min       = -pi;
  theta_candidate = 0;
  numpts          = 16;
  delta_theta     = (theta_max - theta_min) / numpts;
  iter = 0;
  subiter = 0;
  A748 = Algo748();
  while (delta_theta > tol && iter < max_iter)
    thetas = linspace(theta_min, theta_max, numpts);
    residuals = arrayfun(func, thetas);
    [~, min_idx] = min(residuals);
    theta_candidate = thetas(min_idx);
    theta_min       = theta_candidate-delta_theta;
    theta_max       = theta_candidate+delta_theta;
    delta_theta     = (theta_max - theta_min) / numpts;
    iter = iter + numpts;
    numpts = max(4,numpts/2);
    type_min   = ftype(theta_min);
    type_max   = ftype(theta_max);
    type       = ftype(theta_candidate);
    d_func_min = d_func(theta_min);
    d_func_max = d_func(theta_max);
    if (type_min(1) == type_max(1) && type_min(2) == type_max(2) && (d_func_min) * (d_func_max) <= 0) && type_min(1) == type(1) && type_min(2) == type(2)
      [root, subiter] = A748.eval(theta_min, theta_max, d_func);
      tot_iter = iter + subiter;
      %fprintf("Mixed Strategy (pattern + Algo748): Solution found at iter %d\n",tot_iter);
      return;
    end
  end
  root = theta_candidate;
  tot_iter = iter + subiter;
  %fprintf("Mixed Strategy (pattern): Solution found at iter %d\n",tot_iter);
end

%%
function [root, tot_iter] = mixed_strategy_clust2(func, d_func, ftype, tol, max_iter)
  % Initialize variables
  theta_max   = +pi + 2*pi/16;
  theta_min   = -pi - 2*pi/16;
  numpts      = 18; %16
  thetas      = linspace(theta_min, theta_max, numpts);
  Ls          = arrayfun(func, thetas);
  DLs         = arrayfun(d_func, thetas);
  Types       = arrayfun(ftype, thetas, 'UniformOutput', false);
  iter        = numpts;
  clusters    = {};
  [F_theta_candidate, min_idx] = min(Ls);
  theta_candidate = thetas(min_idx);
  % check the pattern high low high of the Ls
  for i = 2:length(Ls)-1
    im = i-1;
    ip = i+1;
    % if i==1
    %   im = length(Ls);
    %   if Ls(im) > Ls(i) && Ls(i) < Ls(ip)
    %     fprintf("test\n")
    %   end
    % end
    if (Ls(im) >= Ls(i) && Ls(i) <= Ls(ip))
      % fprintf("Cerca a i = %d\n",i);
      tmp_struct.theta = [thetas(im), thetas(i), thetas(ip)];
      tmp_struct.L     = [Ls(im), Ls(i), Ls(ip)];
      tmp_struct.dL    = [DLs(im), DLs(i), DLs(ip)];
      tmp_struct.type  = [Types(im), Types(i), Types(ip)];
      % check if max of cluster is above maximum of another eliminate
      % if ~isempty(clusters) && min(tmp_struct.L) > max(cellfun(@(c) max(c.L), clusters))
      %   continue;
      % end
      clusters{end+1} = tmp_struct;
    end
  end
  %
  tot_iter = iter;
  %

  F_root = F_theta_candidate;


  queue = clusters; % Initialize the queue with the initial clusters

  while ~isempty(queue)
    cluster = queue{1};
    queue(1) = []; % Dequeue the first cluster

    if cluster.theta(1) > pi && cluster.theta(2) > pi && cluster.theta(3) > pi
      cluster.theta = cluster.theta-2*pi;
    end
    if cluster.theta(1) < -pi && cluster.theta(2) < -pi && cluster.theta(3) < -pi
      cluster.theta = cluster.theta+2*pi;
    end

    % Check if the cluster has the same type of solution and the derivative at the border have opposite signs
    if cluster.type{1}(1) == cluster.type{2}(1) && cluster.type{2}(1) == cluster.type{3}(1) && ...
       cluster.type{1}(2) == cluster.type{2}(2) && cluster.type{2}(2) == cluster.type{3}(2) && ...
       (cluster.dL(1)) * (cluster.dL(3)) <= 0
      A748 = Algo748();
      [root, subiter] = A748.eval(cluster.theta(1), cluster.theta(3), d_func);
      if root > pi
        root = root-2*pi;
      elseif root < -pi
        root = root+2*pi;
      end
      tot_iter = tot_iter + subiter;
      F_root = func(root);
      if F_root < F_theta_candidate
        theta_candidate = root;
        F_theta_candidate = F_root;
      end
      continue;
    end

    % From cluster compute additional points between theta(1) and theta(3)
    thetas = [cluster.theta(1), (cluster.theta(1) + cluster.theta(2)) / 2, cluster.theta(2), (cluster.theta(2) + cluster.theta(3)) / 2, cluster.theta(3)];
    Ls     = [cluster.L(1), func(thetas(2)), cluster.L(2), func(thetas(4)), cluster.L(3)];
    DLs    = [cluster.dL(1), d_func(thetas(2)), cluster.dL(2), d_func(thetas(4)), cluster.dL(3)];
    Types  = {cluster.type{1}, ftype(thetas(2)), cluster.type{2}, ftype(thetas(4)), cluster.type{3}};
    % [F_theta_candidate, min_idx] = min(Ls);
    % theta_candidate = thetas(min_idx);
    iter = 2;

    new_clusters = {};
    for i = 2:length(Ls)-1
      if Ls(i-1) >= Ls(i) && Ls(i) <= Ls(i+1)
        tmp_struct.theta = [thetas(i-1), thetas(i), thetas(i+1)];
        tmp_struct.L     = [Ls(i-1), Ls(i), Ls(i+1)];
        tmp_struct.dL    = [DLs(i-1), DLs(i), DLs(i+1)];
        tmp_struct.type  = [Types(i-1), Types(i), Types(i+1)];
        new_clusters{end+1} = tmp_struct;
      end
    end

    tot_iter = tot_iter + iter;

    for i = 1:length(new_clusters)
      i_cluster = new_clusters{i};
      if (abs(i_cluster.theta(1) - i_cluster.theta(2)) < tol || abs(i_cluster.theta(2) - i_cluster.theta(3)) < tol) || tot_iter > max_iter
        F_root = i_cluster.L(2);
        root = i_cluster.theta(2);
        subiter = 0;
        if F_root < F_theta_candidate
          theta_candidate = root;
          F_theta_candidate = F_root;
        end
      else
        queue{end+1} = i_cluster; % Enqueue the new cluster
        subiter = 0; % No iteration in this loop step
      end
      tot_iter = tot_iter + subiter;
    end
  end

  root = theta_candidate;
end


%%

function res = sim(mesh, id, figure_flag, save_flag)

  d      = 3;
  x0     = -d;
  y0     = 0;
  xM     = 0;
  yM     = 0;
  xf     = d;
  yf     = 0;
  theta0    = mesh(1,id);
  thetaf    = mesh(2,id);
  k_max     = mesh(3,id);
  theta_int = mesh(4,id);
  xM = 0.9*d*cos(theta_int);
  yM = 0.9*d*sin(theta_int);


  L  = @(thetaM) len_Dubins(x0, y0, theta0, xM, yM, thetaM, k_max ) + ...
                len_Dubins(xM, yM, thetaM, xf, yf, thetaf, k_max );
  DL = @(thetaM) len_Dubins_DR(x0, y0, theta0, xM, yM, thetaM, k_max ) + ...
                len_Dubins_DL(xM, yM, thetaM, xf, yf, thetaf, k_max );

  typeS = @(thetaM) [ type_Dubins(x0, y0, theta0, xM, yM, thetaM, k_max ) ; ...
                      type_Dubins(xM, yM, thetaM, xf, yf, thetaf, k_max ) ];

  tic();
  [thetaM0,n_iter] = mixed_strategy_clust2(L, DL, typeS, 1e-6, 10000);
  elapsed = toc();

  if mod(id,10)==0
    fprintf('------------------------------------------------------------\n');
    fprintf('Test %d/%d\n',id,length(mesh));
    fprintf('  %.2f %% completed\n',100*id/length(mesh));
    fprintf('  theta0 = %g, thetaf = %g, k_max = %g\n',theta0,thetaf,k_max);
    fprintf("  Elapsed time for pattern search: %g\n",elapsed);
    fprintf('  Solution found in N_iter = %d is:  thetaM = %g, L = %g, dL = %g\n',n_iter, thetaM0,L(thetaM0),DL(thetaM0));
    fprintf('------------------------------------------------------------\n');
  end

  theta  = thetaM0;
  n      = n_iter;
  ELLE   = L(thetaM0);
  DELLE  = DL(thetaM0);
  tmp    = typeS(thetaM0);
  TA     = tmp(1);
  TB     = tmp(2);

  res = [theta,n,ELLE,DELLE,TA,TB,elapsed];

  if figure_flag
    clf(gcf);
    DB_A = Dubins();
    DB_B = Dubins();

    DB_A.build( x0, y0, theta0,  xM, yM, thetaM0, k_max );
    DB_B.build( xM, yM, thetaM0, xf, yf, thetaf,  k_max );

    subplot(3,1,1);
    title(gca,sprintf('Dubins75: %d, alpha = %4.3f, beta = %4.3f, kappa = %4.3f',id,theta0,thetaf,k_max));
    hold on;
    DB_A.plot();
    DB_B.plot();
    plot([x0,xM,xf],[y0,yM,yf],'o','MarkerSize',15,'MarkerFaceColor','red');
    hold on
    theta_int_p = linspace(0,2*pi,100);
    xM_p = 0.9*d*cos(theta_int_p);
    yM_p = 0.9*d*sin(theta_int_p);
    plot(xM_p,yM_p,'r--', 'LineWidth',1)
    axis equal
    grid on

    npts     = 1000;
    epsilon  = 1e-4;
    thetas   = linspace(-pi,+pi,npts);
    LAB      = arrayfun(L,  thetas);
    D_LAB    = arrayfun(DL, thetas);
    D_LAB_FD = (arrayfun(L, thetas+epsilon) - arrayfun(L, thetas-epsilon))./(2*epsilon);
    IDX = find( abs(D_LAB_FD) > 5 );
    D_LAB_FD(IDX) = NaN;


    subplot(3,1,2);
    plot(thetas,LAB,'LineWidth',3,'HandleVisibility','off');
    hold on
    plot(thetaM0,L(thetaM0),'o','MarkerSize',15,'MarkerFaceColor','blue');
    [min_tmp, min_idx] = min(LAB);
    plot(thetas(min_idx),min_tmp,'o','MarkerSize',10,'MarkerFaceColor','green');
    grid on
    legend({'$\min_{pattern}$','$\min_{sample}$'})



    subplot(3,1,3);
    plot(thetas,D_LAB,'LineWidth',3);
    hold on
    plot(thetas,D_LAB_FD,'LineWidth',2);
    plot(thetaM0,DL(thetaM0),'o','MarkerSize',10,'MarkerFaceColor','green');
    grid on
    legend({'$DL$','$DL_{FD}$','$\min$'})

    % save the figure
    % create a figure name progressively  i.e _001, _002, ... _999
    %

    if save_flag
      figname = sprintf('./img/Clothoids_test_Dubins75_%03d',id);
      saveas(gcf,figname,'png');
    else
      pause(0.01)
    end

  end




end


