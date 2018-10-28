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
theta0 = 0*1.433674+pi;

S = ClothoidList();
S.push_back( x0, y0, theta0, 0, 0, 608.632424 );
S.push_back( -0.073376, 0, 21.665602 ); 
S.push_back(         0, 0, 13.967496 ); 
S.push_back(  0.072339, 0, 26.803483 ); 
S.push_back(         0, 0, 76.833064 );
S.push_back( -0.003161, 0, 106.008819 ); 
S.push_back(         0, 0, 108.027343 );
S.push_back( -0.003379, 0, 256.145859 );
S.push_back(         0, 0, 7.323885 );
S.push_back( -0.002527, 0, 216.435583 ); 
S.push_back(         0, 0, 233.727813 ); 
S.push_back( -0.000884, 0, 43.695944 ); 
S.push_back(         0, 0, 97.237247 ); 
S.push_back(  0.045987, 0, 26.017783 );
S.push_back(         0, 0, 18.662756 );
S.push_back( -0.049279, 0, 22.635693 ); 
S.push_back(         0, 0, 3.884634 ); 
S.push_back(  0.002349, 0, 94.618213 ); 
S.push_back(         0, 0, 15.677035 );
S.push_back(  0.000723, 0, 44.030556 ); 
S.push_back(         0, 0, 48.016891 ); 
S.push_back( -0.001809, 0, 32.895447 ); 
S.push_back(         0, 0, 50.501013 ); 
S.push_back( -0.013539, 0, 132.829413 ); 
S.push_back(         0, 0, 234.043754 ); 
S.push_back( -0.024510, 0, 46.968166 ); 
S.push_back(         0, 0, 174.535602 );
S.push_back(  0.000057, 0, 93.480793 ); 
S.push_back(         0, 0, 74.093935 ); 
S.push_back(  0.002246, 0, 107.246465 ); 
S.push_back(         0, 0, 24.460791 ); 
S.push_back( -0.000535, 0, 109.331749 );
S.push_back(         0, 0, 149.577842 ); 
S.push_back(  0.000143, 0, 200.204977 );
S.push_back(         0, 0, 87.838354 );
S.push_back(  0.022437, 0, 42.688452 );
S.push_back(         0, 0, 30.118544 ); 
S.push_back( -0.009662, 0, 109.031756 );
S.push_back(         0, 0, 9.262443 ); 
S.push_back(  0.014180, 0, 58.999275 ); 
S.push_back(         0, 0, 924.063894 ); 
S.push_back( -0.011584, 0, 116.757852 ); 
S.push_back(         0, 0, 3.284369 );
S.push_back( -0.013164, 0, 36.956285 ); 
S.push_back(         0, 0, 4.281928 );
S.push_back( -0.005815, 0, 204.111573 ); 
S.push_back(         0, 0, 172.966780 );
S.push_back( -0.000931, 0, 123.006586 );
S.push_back_G1( x0, y0, theta0 ); 
%S.push_back(         0, 0, 309.601534 ); 
S.plot();
axis equal
