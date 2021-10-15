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
theta0 = pi+pi/4;

S = ClothoidList();

addseg  = @(a,r) S.push_back( -1/r, 0, abs(a*r) );
addseg1 = @(l) S.push_back( 0, 0, l );

% STRAIGHT LINE      START GRID
%% STRAIGHT LINE														      START GRID
S.push_back( x0, y0, theta0, 0, 0, 257 ); % first segment

%  CORNER  (1. Variante del parco, dx)
addseg(1.17560115144983,33.6446932777513);

% STRAIGHT LINE
addseg1(69.3886020301746);

%  CORNER (2 ,cambio direzione sx)
addseg(1.08225183870547,-28.0372443981261);

% STRAIGHT LINE
addseg1(156.891281546505);

% CORNER (3. Curva del Rio , dx)
addseg(1.12076347287888,56.0744887962522);

% STRAIGHT LINE
addseg1(177.150952012352);

% CORNER (4 ,cambio direzione dx )
addseg(2.14495670467292,16.8223466388757);

% STRAIGHT LINE
addseg1(54.2676483865244);

% CORNER (5 ,cambio direzione dx ) ...
addseg(1.28777721372631,22.4297955185009);

% STRAIGHT LINE
addseg1(66.2725918538787);

% CORNER (6 ,cambio direzione sx)
addseg(1.19008632273402,-39.2521421573765);

% STRAIGHT LINE
addseg1(334.525723771273);

% CORNER: (congiungente sx)
addseg(0.192440673730988,-123.363875351755);

% STRAIGHT LINE
addseg1(255.294862398743);

% CORNER: (7. Quercia , sx)
addseg(2.4598907976498,-35.8876728296014);

% STRAIGHT LINE
addseg1(83.8848390576514);

% CORNER: (congiungente sx)
addseg(0.181868025357458,-100.934079833254);

% STRAIGHT LINE
addseg1(213.852827215533);
%%
% CORNER: (congiungente dx)
addseg(0.497651064634048,44.8595910370017);

% STRAIGHT LINE
addseg1(53.4220652374333);

% CORNER: (8. Tramonto , dx)
addseg(2.61087521037001,33.6446932777513);

% STRAIGHT LINE
addseg1(52.7112234886856);

% CORNER: (congiungente dx)
addseg(0.151521126772546,67.2893865555026);

% STRAIGHT LINE
addseg1(481.340920785153);

% CORNER: (9. Curvone , dx)
addseg(0.558770009523482,106.541528712879);

% STRAIGHT LINE
addseg1(227.492188959158);

% CORNER: (10 ,cambio direzione dx)
addseg(0.559644334132245,95.3266309536287);

% STRAIGHT LINE
addseg1(118.653844587077);

% CORNER: (11 ,cambio direzione dx)
addseg(1.11496041368553,39.2521421573765);

% STRAIGHT LINE
addseg1(123.384733469072);

% CORNER: (12. Curva del Carro , dx)
addseg(2.69586365772195,23.5512852944259);

% STRAIGHT LINE
addseg1(150.56420705133);

% CORNER: (13 ,cambio di direzione sx)
addseg(1.1002346592359,-39.2521421573765);

% STRAIGHT LINE
addseg1(223.404077635714);

% CORNER: (14 ,cambio direzione sx)
addseg(1.42844614134274,-26.915754622201);

% STRAIGHT LINE
%addseg1(283.563770397571);
S.push_back_G1( x0, y0, theta0); % close curve
%%

S.plot();
axis equal
