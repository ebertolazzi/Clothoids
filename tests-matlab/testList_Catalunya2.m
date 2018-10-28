addpath('../matlab');

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

% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||   START GRID   |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||	      
S.push_back( x0, y0, theta0, 0, 0, 683.18421723 ); % 1
%addseg(683.18421723,0);

% CORNER Elf
addseg(66.00756336,1/41.6943338);

% STRAIGHT LINE
addseg(26.29065828,0);

% CORNER Renault
addseg(66.4510856,-1/55.08714403);

% STRAIGHT LINE
addseg(82.79343713,0);

% CORNER 1 
addseg(198.08208576,1/111.8800341);

% CORNER 2 
addseg(89.40404942,1/162.32043025);

% CORNER 3 
addseg(110.57821361,1/251.9645594);

% STRAIGHT LINE
addseg(225.79185484,0);

% CORNER Repsol 1
addseg(67.69707522,1/52.77611157);

% ORNER Repsol 2
addseg(56.6410650,1/68.88744017);

% CORNER Repsol 3
addseg(111.76662754,1/107.07259597);

% STRAIGHT LINE
addseg(179.11072228,0);

% CORNER Seat
addseg(91.12679306,-1/35.58186104);

% STRAIGHT LINE
addseg(113.93721138,0);

% CORNER 1
addseg(53.88581637,-1/469.97031844);

% CORNER 2
addseg(84.29604548,-1/186.41064472);

% STRAIGHT LINE
addseg(93.10471419,0);

% CORNER
addseg(69.50876223,-1/42.15162403);

% STRAIGHT LINE
addseg(23.14519818,0);

% CORNER
addseg(49.5041283,1/74.87247217);

% STRAIGHT LINE
addseg(176.26004018,0);

% CORNER Campsa
addseg(143.1226559,1/94.70681912);

% STRAIGHT LINE
addseg(517.7104027,0);

% CORNER la caixa 1  
addseg(52.34464597,-1/35.52471376);

% CORNER la caixa 2  
addseg(32.30604152,-1/51.12970729);

% CORNER la caixa 3  
addseg(86.85492868,-1/105.23505766);

% CORNER la caixa 4  
addseg(53.41214976,-1/120.83961628);

% STRAIGHT LINE
addseg(55.07426982,0);

% CORNER Banc Sabadel
addseg(184.02689381,1/58.30616547);

% STRAIGHT LINE
addseg(144.69506316,0);

% CORNER new holland 
addseg(96.42612713,1/75.62996779);

% STRAIGHT LINE
addseg(165.30277047,0);

% CORNER 
addseg(137.27937829,1/88.9444796);

% STRAIGHT LINE
%addseg(360.27575009,0);
S.push_back_G1(x0,y0,theta0);

S.plot();
axis equal