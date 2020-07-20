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
theta0 = 0*pi/2;

S = ClothoidList();

to_rad  = @(a) a*pi/180;
addseg  = @(a,r) S.push_back( 1/r, 0, abs(to_rad(a)*r) );
addseg1 = @(l) S.push_back( 0, 0, l );

S.push_back( x0, y0, theta0, 0, 0, 378.364557 ); % first segment

addseg( 44.000000,-339.539241 );
addseg( 15.000000,-60.744304 );
addseg1( 38.324051 );
addseg( 60.000000, 18.808861 );
addseg1( 13.408861 );
addseg( 52.000000, -63.341772 );
addseg1( 517.221519 );
addseg( 80.000000,-29.825316 );
addseg1( 16.906329 );
addseg( 93.000000, -29.916456 );
addseg1( 305.510127 );
addseg( 24.000000,56.734177);
addseg( 76.000000,20.517722);
addseg( 77.000000,58.545570);
addseg1( 28.720253 );
addseg( 18.000000,111.782278 ); %%%% +
addseg1( 223.177215 );
addseg( 80.526316,-28.013924 );
addseg1( 22.910127 );
addseg( 84.000000,-28.913924);
addseg1( 54.341772);
addseg( 19.000000,-74.255696); %%% -
addseg1( 655.821519);
addseg( 86.000000,37.321519);
addseg( 64.000000,-56.039241);
addseg1( 207.945570);
addseg( 63.000000,-44.532911);
addseg( 60.000000,-60.539241);
addseg1( 34.382278);
addseg( 62.000000,20.312658);
addseg1( 10.708861);
addseg( 59.000000,26.213924);
addseg1( 126.786076);
addseg( 88.000000,-22.613924);
addseg1( 41.365823);
addseg( 78.000000,-24.208861);
addseg1( 170.724051);


%addseg1(283.563770397571);
%S.push_back_G1( x0, y0, theta0); % close curve
%%

S.plot();
axis equal