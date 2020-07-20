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

C = ClothoidCurve();
C.build_G1( 0, 0, -pi, 1, 3, pi );

P = PolyLine();
P.approx( C, 0.001 );

C.plot(1000,'LineWidth',3);
hold on;
[ x, y ] = P.polygon();
P.plot();

axis equal

P.delete();
C.delete();
