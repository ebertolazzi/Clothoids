addpath('../matlab');

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

x0 = 0;
y0 = 0;
x1 = 1;
y1 = 0;
kk = 1;

thmin = 0.5*pi;
thmax = 0.8*pi;

r0=thmin:(thmax-thmin)/2:thmax;

thr0={r0,r0};
thr1={r0,-r0};

L = 0.4;
LAB = {'(a)','(b)','(c)','(d)'};

ba = Biarc();

aa = 0.04;
bb = 0.5-2*aa;

figure('Position',[ 1 1 600 600]);


for kk=[0,1]

  switch(1+2*kk)
  case 1; subplot('Position',[aa aa bb bb]);
  case 2; subplot('Position',[aa+0.5 aa bb bb]);
  case 3; subplot('Position',[aa+0.5 aa+0.5 bb bb]);
  case 4; subplot('Position',[aa aa+0.5 bb bb]);
  end

  hold off
  plot([x0,x1],[y0,y1],'k');
  hold on

  for theta1=thr1{kk+1}
    for theta0=thr0{kk+1}
      ok = ba.build(x0,y0,theta0,x1,y1,theta1);
      if ok
        ba.plot();
        quiver(x0,y0,L*cos(theta0),L*sin(theta0),'LineWidth',2);
        quiver(x1,y1,L*cos(theta1),L*sin(theta1),'LineWidth',2);
      end
    end
  end
  
  axis equal
  set(gca,'DataAspectRatio',[1,1,1]);
  set(gca,'XTick',[-2,-1,0,1,2]);
  set(gca,'XTickLabel',{'-2','-1','0','1','2'});
  set(gca,'YTick',[-2,-1,0,1,2]);
  set(gca,'YTickLabel',{'-2','-1','0','1','2'});
  title([LAB{1+2*kk} ' Present Method'])

  switch(2+2*kk)
  case 1; subplot('Position',[aa aa bb bb]);
  case 2; subplot('Position',[aa+0.5 aa bb bb]);
  case 3; subplot('Position',[aa+0.5 aa+0.5 bb bb]);
  case 4; subplot('Position',[aa aa+0.5 bb bb]);
  end

  hold off
  plot([x0,x1],[y0,y1],'k');
  hold on

  for theta1=thr1{kk+1}
    for theta0=thr0{kk+1}
      a0 = theta0+pi/2;
      a1 = theta1+pi/2;
      p=[x0,x1;y0,y1];
      u=[cos(a0),cos(a1);sin(a0),sin(a1)];
      bi_arc = fnrfn(rscvn(p,u),[0.5,1.5]);
      bspline_plot(bi_arc,false);
      ok = ba.build(x0,y0,theta0,x1,y1,theta1);
      if ok
        ba.plot();
        quiver(x0,y0,L*cos(theta0),L*sin(theta0),'LineWidth',2);
        quiver(x1,y1,L*cos(theta1),L*sin(theta1),'LineWidth',2);
      end
    end
  end

  axis equal
  set(gca,'DataAspectRatio',[1,1,1]);
  set(gca,'XTick',[-2,-1,0,1,2]);
  set(gca,'XTickLabel',{'-2','-1','0','1','2'});
  set(gca,'YTick',[-2,-1,0,1,2]);
  set(gca,'YTickLabel',{'-2','-1','0','1','2'});
  title([LAB{2+2*kk} ' Matlab rscvn'])
  
end

if false
  matlab2tikz('figure3.tex', ...
              'standalone',true, ...
              'extraaxisoptions',{'xlabel style={font=\LARGE}','ylabel style={font=\LARGE}','ticklabel style={font=\LARGE}'}, ...
              'extraTikzpictureOptions',{'cap=round'}); 
end

