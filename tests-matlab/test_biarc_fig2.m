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
choose = 1;

thmin = 0;
thmax = 0.5*pi;

LAB = {'(a)','(b)','(c)','(d)'};

ba = Biarc();

aa = 0.04;
bb = 0.5-2*aa;

figure('Position',[ 1 1 600 900]);

for theta1=thmin:(thmax-thmin)/3:thmax

  switch(choose)
  case 1; subplot('Position',[aa aa bb bb]);
  case 2; subplot('Position',[aa+0.5 aa bb bb]);
  case 3; subplot('Position',[aa+0.5 aa+0.5 bb bb]);
  case 4; subplot('Position',[aa aa+0.5 bb bb]);
  end

  hold off
  plot([x0,x1],[y0,y1],'k');
  hold on

  for theta0=[-pi:pi/20:pi]
    ok = ba.build(x0,y0,theta0,x1,y1,theta1);
    if ok
      ba.plot();
    end
  end
  axis([-1,1.2,-1.4,1.4]);
  axis equal
  xlabel(LAB{choose});

  set(gca,'DataAspectRatio',[1,1,1]);
  set(gca,'XTick',[-2,-1,0,1,2]);
  set(gca,'XTickLabel',{'-2','-1','0','1','2'});
  set(gca,'YTick',[-2,-1,0,1,2]);
  set(gca,'YTickLabel',{'-2','-1','0','1','2'});
  legend({},'FontSize',12);
  legend('off');
  switch(choose)
  case 1; title('\vartheta_1 = 0');
  case 2; title('\vartheta_1 = \pi/6');
  case 3; title('\vartheta_1 = \pi/3');
  case 4; title('\vartheta_1 = 2\pi/3');
  end
  choose = choose+1;
end

if false
  matlab2tikz('figure2.tex', ...
              'standalone',true, ...
              'extraaxisoptions',{'xlabel style={font=\LARGE}','ylabel style={font=\LARGE}','ticklabel style={font=\LARGE}'}, ...
              'extraTikzpictureOptions',{'cap=round'}); 
end

ba.delete();
