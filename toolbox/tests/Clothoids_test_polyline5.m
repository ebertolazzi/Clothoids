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
     2.953123242,0.925003242,-0.915616758,-2.989063158,-3.178123242,...
     -3.178123242, -1.707808758 ];


SG2 = ClothoidSplineG2();
S   = SG2.buildP2( X, Y );

S1 = S.copy();
S1.rotate(pi,0,0);
S1.translate(0,-2);

S.rotate(pi/2,0,0);
S1.rotate(pi/2,0,0);

P1 = PolyLine();
P  = PolyLine();

P1.approx( S1, 0.1 );
P.approx( S, 0.1 );

if true

P1.plot();
hold on;
P.plot();

tic
[ s1, s2 ] = P.intersect( P1 );
toc

XY1 = P.eval( s2 );
XY2 = P1.eval( s1 );

if ~isempty(s1)
  plot( XY1(1,:), XY1(2,:), 'ob', 'LineWidth', 2, 'Markersize', 10 );
end

if ~isempty(s2)
  plot( XY2(1,:), XY2(2,:), 'oc', 'LineWidth', 2, 'Markersize', 8 );
end

else  
    
P1.plot();
hold on;
S.plot();

tic
[ s3, s4 ] = S.intersect( P1 );
toc

XY3 = S.eval( s3 );
XY4 = P1.eval( s4 );

if ~isempty(s3)
  plot( XY3(1,:), XY3(2,:), 'ob', 'LineWidth', 2, 'Markersize', 10 );
end

if ~isempty(s4)
  plot( XY4(1,:), XY4(2,:), 'oc', 'LineWidth', 2, 'Markersize', 8 );
end

end

axis equal
