%=========================================================================%
%                                                                         %
%  Autor: Enrico Bertolazzi                                               %
%         Department of Industrial Engineering                            %
%         University of Trento                                            %
%         enrico.bertolazzi@unitn.it                                      %
%                                                                         %
%=========================================================================%
% Driver test program to check Clothoids lib                              %
%=========================================================================%

close all;

X = [2.9265642,2.6734362,2.5109322,1.9078122,1.1859282,1.9249962, ...
     2.8265562,0.00468420000000025,-2.826567,-1.9437558,-1.1859438, ...
     -1.9062558,-2.501565,-2.6734386,-2.9265642,-2.6187522,-1.1406318, ...
     -0.8968758,-1.4562558,-1.9062558,-0.00468780000000013,1.9078122, ...
     1.4468682,0.8968722,1.1406282,2.6187522, 2.9265642 ];
Y = [-1.707808758,-1.707808758,-2.367185958,-2.582810358,-2.582810358, ...
     -1.167184758,0.915619242,3.178123242,0.915619242,-1.150000758, ...
     -2.582810358,-2.582810358,-2.393750358,-1.707808758,-1.707808758, ...
     -3.178123242,-3.178123242,-2.989063158,-0.915616758,0.925003242, ...
     2.953123242,0.925003242,-0.915616758,-2.989063158,-3.178123242,-3.178123242, -1.707808758 ];

SC = ClothoidSplineG2();
S1 = ClothoidList();
S2 = ClothoidList();
S3 = ClothoidList();

subplot(3,1,1);

S = SC.buildP2( X, Y );
S1.copy(S);
S2.copy(S);
S3.copy(S);

fmt1 = {'Color','red','LineWidth',4};
fmt2 = {'Color','blue','LineWidth',4};

x0 = 0;
y0 = 0;

subplot(2,2,1);
S.plot( 100, fmt1, fmt2);
[ x1, y1, s1, t1, iflag1, dst1 ] = S.closestPoint( x0, y0 );
%[ x, y, s, t, iflag, dst ] = S.closestPointInSRange( PNT(1), PNT(2), s_begin, s_end );
hold on;
%[x1,y1] = S.eval( s );
plot([x0,x1],[y0,y1],'-s','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','red');
axis equal
title('one');

subplot(2,2,2);
s_min = 0;
s_max = 10;
S1.trim(s_min, s_max);
S1.plot( 100, fmt1, fmt2);
res1 = S1.closestPoint( x0, y0 );
res2 = S.closestPointInSRange( x0, y0, s_min, s_max );
idx  = S.s_to_index( s_min );
hold on;
%[x1,y1] = S1.eval( s );
plot([x0,res1.x],[y0,res1.y],'-s','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','red');
plot([x0,res2.x],[y0,res2.y],'-o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b');
axis equal
fprintf('iflag+idx = %d, icurve = %d\n\n', res1.iflag+idx, res2.icurve );
title('two');

subplot(2,2,3);
s_min = 10;
s_max = 43;
S2.trim( s_min, s_max );
S2.plot( 100, fmt1, fmt2);
res1 = S2.closestPoint( x0, y0 );
res2 = S.closestPointInSRange( x0, y0, s_min, s_max );
idx  = S.s_to_index( s_min );
hold on;
%[x1,y1] = S1.eval( s );
plot([x0,res1.x],[y0,res1.y],'-s','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','red');
plot([x0,res2.x],[y0,res2.y],'-o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b');
axis equal
fprintf('iflag+idx = %d, icurve = %d\n\n', res1.iflag+idx, res2.icurve );
title('three');

subplot(2,2,4);
s_min = 35;
s_max = 5;
S3.trim(s_min, s_max);
S3.plot( 100, fmt1, fmt2);
res1 = S3.closestPoint( x0, y0 );
res2 = S.closestPointInSRange( x0, y0, s_min, s_max );
idx  = S.s_to_index( s_min );
hold on;
%[x1,y1] = S1.eval( s );
plot([x0,res1.x],[y0,res1.y],'-s','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','red');
plot([x0,res2.x],[y0,res2.y],'-o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b');
axis equal
fprintf('iflag+idx = %d, icurve = %d\n\n', res1.iflag+idx, res2.icurve );
title('four');

