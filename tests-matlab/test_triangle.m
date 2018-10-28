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

addpath('../matlab');

% check constructors
x0     = 0;
y0     = 2;
x1     = 4;
y1     = 3.5;
x2     = 6;
y2     = 1;

T1 = Triangle2D( x0, y0, x1, y1, x2, y2 );
T2 = Triangle2D( [x0, y0], [x1, y1], [x2, y2] );
T2.translate(5,6);

T1.plot('red');
hold on;
T2.plot('blue','Color','red','Linewidth',3);

axis equal;
