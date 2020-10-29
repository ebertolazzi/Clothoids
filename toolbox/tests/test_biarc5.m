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

x0 = -1;
y0 = 0;
x1 = 1;
y1 = 0;

hold off
plot([x0,x1],[y0,y1],'k');
hold on

for theta1=[-pi:pi/5:pi]

  for theta0=[-pi:pi/5:pi]
    L1 = Biarc(x0,y0,theta0,x1,y1,theta1);
    L1.plot();
  end
end

xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
axis equal
