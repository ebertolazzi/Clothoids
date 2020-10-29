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
x0     = [-4,-4,-4,-4];
y0     = [6,6,6,6];
x1     = [2,3,4,4];
y1     = [1,-1,2,3];
theta0 = 0;
theta1 = [0,pi/2,-pi/2,pi/4];

aa = 0.04;
bb = 0.5-2*aa;

figure('Position',[ 1 1 800 800]);

for kk=1:4

  switch(kk)
  case 1; subplot('Position',[aa aa bb bb]);
  case 2; subplot('Position',[aa+0.5 aa bb bb]);
  case 3; subplot('Position',[aa+0.5 aa+0.5 bb bb]);
  case 4; subplot('Position',[aa aa+0.5 bb bb]);
  end
  
  S = ClothoidList();
  iter = S.build_3arcG2( x0(kk), y0(kk), theta0, 0, x1(kk), y1(kk), theta1(kk), 0 );

  S.plot(400,{'Color','blue','LineWidth',3},{'Color','red','LineWidth',3});
  
  N = 10;

  x         = (-10 + 20 * rand(N,1)).';
  y         = (-5  + 20 * rand(N,1)).';
  
  [s,t]     = S.find_coord(x,y)

  [xx,yy]   = S.eval(s,t);
  [xxx,yyy] = S.eval(s);
  

  hold on;  
  plot( x,   y, 'ro', 'MarkerSize',20,'MarkerEdgeColor','blue','MarkerFaceColor',[0.5,0.5,0.5]);
  plot( xx, yy, 'b*', 'MarkerSize',15,'MarkerEdgeColor','red','MarkerFaceColor',[0.9,0.9,0.5]);
  plot( [x;xxx;xx], [y;yyy;yy] );
  
  sss = (0:S.length()/1000:S.length()).';
  [x,y] = S.eval( sss, 1 );
  plot( x, y );
  [x,y] = S.eval( sss, -0.5 );
  plot( x, y );

  axis equal;
  %
end
