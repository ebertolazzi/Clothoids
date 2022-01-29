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

addseg  = @(l,c) S.push_back( -c, 0, l );
addseg1 = @(l) S.push_back( 0, 0, l );

S.push_back( x0, y0, theta0, 0, 0, 20 ); % 1

addseg1(333.56399);
% CORNER
addseg( 94.430509,-0.014788752);

% STRAIGHT LINE
addseg1(288.876735);

% CORNER
addseg(76.993736,-0.033496099);

% STRAIGHT LINE
addseg1(111.792956);

% CORNER
addseg(32.317242,-0.025091069);

% STRAIGHT LINE
addseg1(126.901826);

% CORNER
addseg(55.0546292,0.029835693);

% STRAIGHT LINE
addseg1(84.134941);

% CORNER
addseg(113.685868,53.072559^(-1))

% STRAIGHT LINE
addseg1(192.264108);

% CORNER
addseg(73.591598,-34.35873^(-1));

% STRAIGHT LINE														 
addseg1(292.868903);

% CORNER
addseg(26.988876,-61.753526^(-1));

% STRAIGHT LINE														 
addseg1(134.851756);

% CORNER
addseg(77.107342,-63.684985^(-1));

% CORNER
addseg(54.232517,-36.444969^(-1));

% STRAIGHT LINE
addseg1(105.707786);

% CORNER
addseg(68.235168,-85.046849^(-1));

% STRAIGHT LINE														 
addseg1(56.609801);

% CORNER
addseg(56.141288,67.287775^(-1));

% STRAIGHT LINE
addseg1(56.524823);

% CORNER
addseg(94.573034,35.332293^(-1));

% STRAIGHT LINE
addseg1(281.247571);

% CORNER
addseg(42.124813,35.485691^(-1));

% CORNER
addseg(476.610109,-241.713468^(-1));

% CORNER
addseg(59.981871,-31.11497^(-1));

% STRAIGHT LINE
addseg1(497.63868);
%addseg1(20.0);
S.push_back_G1(x0,y0,theta0); 

S.plot();
axis equal
