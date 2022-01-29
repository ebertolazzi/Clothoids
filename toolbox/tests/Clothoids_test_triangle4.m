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

% check constructors

x0 = 0; y0 = -1;
x1 = 1; y1 = 0;
x2 = 0; y2 = 7;
T1 = Triangle2D( x0, y0, x1, y1, x2, y2 );

x0 = 0; y0 = 0;
x1 = 6; y1 = 0;
x2 = 6; y2 = 3;
T2 = Triangle2D( x0, y0, x2, y2, x1, y1 );

fmt1 = {'EdgeColor','red','Linewidth',2};
fmt2 = {'EdgeColor','green','Linewidth',2};
c1   = {'red','FaceAlpha',0.1};
c2   = {'blue','FaceAlpha',0.1};

figure('Position',[ 1 1 800 800]);

T1.translate(-5,0);
T2.translate(-2,0);

for k=1:100

  T1.translate(0.1,0);

  hold off;
  if T1.overlap(T2)
    T1.plot2(c1{:},fmt1{:});
    hold on;
    T2.plot2(c2{:},fmt1{:});
  else
    T1.plot2(c1{:},fmt2{:});
    hold on;
    T2.plot2(c2{:},fmt2{:});
  end

  plot(-10,-10,'o');
  plot(10,10,'o');

  axis equal;
  drawnow;
end