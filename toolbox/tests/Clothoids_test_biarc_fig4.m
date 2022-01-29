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

format long;
addpath(fullfile(pwd,'../matlab2tikz','src'));

close all;

x0=0;
y0=0;
x1=1;
y1=0;

theta0 = pi/2;
theta1 = pi/2;
a0     = theta0+pi/2;
a1     = theta1+pi/2;

%AXE = [-0.02 1.02 -0.52 0.52];
AXE = [0.1,1.1,-0.6,0.6];

L = 0.4;

p=[x0,x1;y0,y1];
u=[cos(a0),cos(a1);sin(a0),sin(a1)];

ba = Biarc();

aa = 0.04;
bb = 0.5-2*aa;

figure('Position',[ 1 1 600 600]);

subplot('Position',[aa+0.5 aa bb bb]);

hold off
plot(p(1,:),p(2,:),'k')
hold on
bi_arc = fnrfn(rscvn(p,u),[0.5,1.5]);
bspline_plot(bi_arc,false);
quiver(x0,y0,L*cos(theta0),L*sin(theta0),'LineWidth',2);
quiver(x1,y1,L*cos(theta1),L*sin(theta1),'LineWidth',2);

axis(AXE);
axis equal
title('(b) Matlab rscvn')

set(gca,'DataAspectRatio',[1,1,1]);
set(gca,'XTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]);
set(gca,'XTickLabel',{'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
set(gca,'YTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]);
set(gca,'YTickLabel',{'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});

subplot('Position',[aa aa bb bb]);

hold off
plot(p(1,:),p(2,:),'k')
hold on
ok = ba.build(x0,y0,theta0,x1,y1,theta1);
if ok
  ba.plot();
  quiver(x0,y0,L*cos(theta0),L*sin(theta0),'LineWidth',2);
  quiver(x1,y1,L*cos(theta1),L*sin(theta1),'LineWidth',2);
end
axis(AXE);
axis equal
title('(a) Present Method');


set(gca,'DataAspectRatio',[1,1,1]);
set(gca,'XTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]);
set(gca,'XTickLabel',{'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
set(gca,'YTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]);
set(gca,'YTickLabel',{'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});

theta0 = pi/2-10000*eps;
theta1 = pi/2-10000*eps;
a0     = theta0+pi/2;
a1     = theta1+pi/2;

p=[x0,x1;y0,y1];
u=[cos(a0),cos(a1);sin(a0),sin(a1)];

subplot('Position',[aa+0.5 aa+0.5 bb bb]);

hold off
plot(p(1,:),p(2,:),'k')
hold on
bi_arc = fnrfn(rscvn(p,u),[0.5,1.5]);
bspline_plot(bi_arc,false);
quiver(x0,y0,L*cos(theta0),L*sin(theta0),'LineWidth',2);
quiver(x1,y1,L*cos(theta1),L*sin(theta1),'LineWidth',2);

axis(AXE);
axis equal
title('(d) Matlab rscvn')

set(gca,'DataAspectRatio',[1,1,1]);
set(gca,'XTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]);
set(gca,'XTickLabel',{'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
set(gca,'YTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]);
set(gca,'YTickLabel',{'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});

subplot('Position',[aa aa+0.5 bb bb]);

hold off
plot(p(1,:),p(2,:),'k')
hold on
ok = ba.build(x0,y0,theta0,x1,y1,theta1);
if ok
  ba.plot();
  quiver(x0,y0,L*cos(theta0),L*sin(theta0),'LineWidth',2);
  quiver(x1,y1,L*cos(theta1),L*sin(theta1),'LineWidth',2);
end
axis(AXE);
axis equal

set(gca,'DataAspectRatio',[1,1,1]);
set(gca,'XTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]);
set(gca,'XTickLabel',{'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
set(gca,'YTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]);
set(gca,'YTickLabel',{'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
title('(c) Present Method');

if false
  matlab2tikz('figure4.tex', ...
              'standalone',true, ...
              'extraaxisoptions',{'xlabel style={font=\LARGE}','ylabel style={font=\LARGE}','ticklabel style={font=\LARGE}'}, ...
              'extraTikzpictureOptions',{'cap=round'}); 
end


