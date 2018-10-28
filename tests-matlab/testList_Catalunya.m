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

% |||||||||||||||||||||||||||||||||||||   START GRID   |||||||||||||||||||||||||||||||||||||||||
S.push_back( x0, y0, theta0, 0, 0, 684.8194564 );
%addseg1(684.8194564);

% CORNER Elf
addseg(61.37799034,1/37.86990612);

% STRAIGHT LINE
addseg1(18.75920636);

% CORNER Renault
addseg(80.069112926,-1/64.36843302);

% STRAIGHT LINE
addseg1(80.68297222);

% CORNER 1 
addseg(63.11609498,1/104.64002904);

% CORNER 2 
addseg(163.8592568,1/115.87061646);

% CORNER 3 
addseg(156.52619938,1/212.94807617);

% STRAIGHT LINE
addseg1(231.8817849);

% CORNER Repsol 1
addseg(105.68477598,1/56.81280378);

% CORNER Repsol 2
addseg(129.95040312,1/100.24537215);

% STRAIGHT LINE
addseg1(180.49228976);

% CORNER Seat
addseg(89.44476175,-1/35.11755175);

% STRAIGHT LINE
addseg1(107.04294776);

% CORNER 1
addseg(60.48237024,-1/469.97031844);

% CORNER 2
addseg(84.29604548,-1/186.41064472);

% STRAIGHT LINE
addseg1(91.23648198);

% CORNER
addseg(68.0044283,-1/42.61284759);

% STRAIGHT LINE
addseg1(28.82157832);

% CORNER
addseg(45.5246094,1/74.87247217);

% STRAIGHT LINE
addseg1(176.26004018);

% CORNER Campsa
addseg(143.1226559,1/94.70681912);

% STRAIGHT LINE
addseg1(515.31737357);

% CORNER la caixa 1  
addseg(59.80535775,-1/35.58231261);

% CORNER la caixa 2  
addseg(51.93254416,-1/69.93813673);

% CORNER la caixa 3  
addseg(113.90487362,-1/119.98115513);

% STRAIGHT LINE
addseg1(50.34381636);

% CORNER Banc Sabadel  
addseg(186.01599897,1/58.9644632);

% STRAIGHT LINE
addseg1(137.73626904);

% CORNER new holland 
addseg(106.47043326,1/83.40969586);

% STRAIGHT LINE
addseg1(149.97055459);

% CORNER 
addseg(155.00375877,1/100.42825682);

% STRAIGHT LINE
%addseg1(349.10204373);
S.push_back_G1( x0, y0, theta0); % close curve

S.plot();
axis equal