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

% check constructors
offs = 8+0.25;
PNT = [     ...
  0,     0; ...
  1,     0; ...
  2,     0; ...
  3,     0; ...
  4,     0; ...
  5,     0; ...
  6,     0; ...
  7,     0; ...
  8,     0; ...
  0+offs,     1; ...
  1+offs,     1; ...
  2+offs,     1; ...
  3+offs,     1; ...
  4+offs,     1; ...
  5+offs,     1; ...
  6+offs,     1; ...
  7+offs,     1; ...
];

theta0 = 0;
theta1 = 0;

figure('Position',[ 1 1 800 800]);

S1 = ClothoidList();
S2 = ClothoidList();
S3 = ClothoidList();
S4 = ClothoidList();

x = PNT(:,1);
y = PNT(:,2);

wmin = -0.1*ones(size(x));
wmax =  0.1*ones(size(x));

ok = S1.build_G2( x, y, theta0, 0, theta1, 0 );
ok = S2.build_G2_with_target( x, y, wmin, wmax, theta0, theta1, 'lenght' );
ok = S3.build_G2_with_target( x, y, wmin, wmax, theta0, theta1, 'curvature' );
ok = S4.build_G2_with_target( x, y, wmin, wmax, theta0, theta1, 'jerk' );

subplot( 4, 1, 1);
len = S1.length();
ss  = 0:len/1000:len;
[xx,yy,theta] = S1.evaluate( ss );

hold on;
plot( x,  y,  'ro', 'MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor',[0.5,0.5,0.5]);
plot( xx, yy, 'b-', 'LineWidth', 2 );
axis equal;

title( 'G2' );

subplot( 4, 1, 2);

len = S2.length();
ss  = 0:len/1000:len;
[xx,yy,theta] = S2.evaluate( ss );
[xp, yp]      = S2.get_XY();

hold on;
plot( xx, yy, 'g-', 'LineWidth', 2 );
plot( x,  y,  'ro', 'MarkerSize',8,'MarkerEdgeColor','blue','MarkerFaceColor',[0.5,0.5,0.5]);
plot( xp, yp, 'yo', 'MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor',[0.5,0.5,0.5]);
axis equal;

title( 'min length' );

subplot( 4, 1, 3);

len = S3.length();
ss  = 0:len/1000:len;
[xx,yy,theta] = S3.evaluate( ss );
[xp, yp]      = S3.get_XY();

hold on;
plot( xx, yy, 'r-', 'LineWidth', 2 );
plot( x,  y,  'ro', 'MarkerSize',8,'MarkerEdgeColor','blue','MarkerFaceColor',[0.5,0.5,0.5]);
plot( xp, yp, 'yo', 'MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor',[0.5,0.5,0.5]);
axis equal;

title( 'min curvature' );

subplot( 4, 1, 4);

len = S4.length();
ss  = 0:len/1000:len;
[xx,yy,theta] = S4.evaluate( ss );
[xp, yp]      = S4.get_XY();

hold on;
plot( xx, yy, 'c-', 'LineWidth', 2 );
plot( x,  y,  'ro', 'MarkerSize',8,'MarkerEdgeColor','blue','MarkerFaceColor',[0.5,0.5,0.5]);
plot( xp, yp, 'yo', 'MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor',[0.5,0.5,0.5]);
axis equal;

title( 'min jerk' );
