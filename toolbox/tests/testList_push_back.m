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


x0     = 0;
y0     = 0;
theta0 = 0;

S = ClothoidList();

addseg  = @(l,c) S.push_back( c, 0, l );
addseg1 = @(l) S.push_back( 0, 0, l );

S.push_back( x0, y0, theta0, 0, 0, 366.933546 ); % 1
addseg( 36.395335, 0.005014);  % 2
addseg1( 17.508761);           % 3
addseg( 26.719483, 0.058793);  % 4
addseg( 62.995808, 0.018212);  % 5
addseg1( 392.057882);          % 6
addseg( 27.236252, -0.010979); % 7
addseg1( 51.298757);           % 8

% check constructors
x0     = 0;
y0     = 2;
theta0 = 0;
k0     = 1/30;
dk     = 0;
L      = 100;

L1 = ClothoidCurve( x0, y0, theta0, k0, dk, L );

S2 = ClothoidList();

S2.info();

subplot(2,1,1)
S2.push_back(S);
S2.plot();

subplot(2,1,2)
S.changeOrigin(S2.xEnd(),S2.yEnd());
S2.push_back(S);
S2.plot();


axis equal
