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
theta0 = pi/3; % 0*1.628490;

S = ClothoidList();
S.push_back( x0, y0, theta0, 0, 0, 41.579688 );
S.push_back( -0.003783, 0, 63.562567 );
S.push_back(         0, 0, 8.305660 ); 
S.push_back( -0.006917, 0, 27.206744 ); 
S.push_back(         0, 0, 48.881450 ); 
S.push_back( -0.047178, 0, 12.849138 ); 
S.push_back(         0, 0, 11.444336 ); 
S.push_back( -0.014110, 0, 29.642836 ); 
S.push_back(         0, 0, 30.892467 ); 
S.push_back(  0.000988, 0, 99.462866 ); 
S.push_back(         0, 0, 24.210783 ); 
S.push_back( -0.004968, 0, 35.347462 ); 
S.push_back(         0, 0, 35.348945 ); 
S.push_back(  0.015191, 0, 21.054050 );
S.push_back(         0, 0, 60.781081 ); 
S.push_back( -0.006323, 0, 44.006296 ); 
S.push_back(         0, 0, 56.399726 ); 
S.push_back(  0.010425, 0, 41.994771 ); 
S.push_back(         0, 0, 13.783875 ); 
S.push_back(  0.018343, 0, 69.965838 );
S.push_back(         0, 0, 6.945343 ); 
S.push_back(  0.013344, 0, 40.409628 ); 
S.push_back(         0, 0, 27.950105 ); 
S.push_back( -0.033342, 0, 45.333484 ); 
S.push_back(         0, 0, 202.269808 ); 
S.push_back( -0.092205, 0, 11.869725 ); 
S.push_back(         0, 0, 2.932659 ); 
S.push_back( -0.041730, 0, 23.952860 ); 
S.push_back(         0, 0, 2.995485 ); 
S.push_back( -0.040513, 0, 9.528251 ); 
S.push_back(         0, 0, 9.006347 ); 
S.push_back(  0.044142, 0, 13.093110 );
S.push_back(         0, 0, 28.103240 ); 
S.push_back( -0.009651, 0, 15.899273 ); 
S.push_back(         0, 0, 5.529594 ); 
S.push_back(  0.094526, 0, 17.664081 ); 
S.push_back(         0, 0, 1.429494 ); 
S.push_back(  0.165028, 0, 10.890622 ); 
S.push_back(         0, 0, 3.754570 ); 
S.push_back( -0.006090, 0, 19.310597 ); 
S.push_back(         0, 0, 28.439009 ); 
S.push_back( -0.068378, 0, 12.774536 ); 
S.push_back(         0, 0, 2.867507 ); 
S.push_back( -0.072714, 0, 13.395234 ); 
S.push_back(         0, 0, 45.671900 ); 
S.push_back( -0.029722, 0, 21.904454 ); 
S.push_back(         0, 0, 10.581918 ); 
S.push_back( -0.153921, 0, 8.718188 ); 
S.push_back(         0, 0, 19.818733 ); 
S.push_back( -0.002107, 0, 113.808848 ); 
S.push_back(         0, 0, 21.238633 ); 
S.push_back( -0.003151, 0, 126.763702 ); 
S.push_back(         0, 0, 14.122690 ); 
S.push_back( -0.002673, 0, 133.328651 ); 
S.push_back(         0, 0, 27.676971 ); 
S.push_back( -0.003514, 0, 101.102754 );
S.push_back(         0, 0, 73.721886 ); 
S.push_back(  0.111113, 0, 10.272486 ); 
S.push_back(         0, 0, 1.169014 );
S.push_back( -0.150817, 0, 7.689617 ); 
S.push_back(         0, 0, 2.250445 ); 
S.push_back( -0.114109, 0, 8.804205 ); 
S.push_back(         0, 0, 1.819195 ); 
S.push_back(  0.120902, 0, 8.527511 ); 
S.push_back(         0, 0, 228.748326 ); 
S.push_back(  0.036398, 0, 22.098581 ); 
S.push_back(         0, 0, 14.805693 ); 
S.push_back(  0.006796, 0, 100.206139 ); 
S.push_back(         0, 0, 21.625315 ); 
S.push_back(  0.023719, 0, 32.517854 ); 
S.push_back(         0, 0, 7.108205 ); 
S.push_back( -0.021477, 0, 31.247850 );
S.push_back(         0, 0, 106.260707 );
S.push_back( -0.058992, 0, 21.978838 );
S.push_back(         0, 0, 7.876341 ); 
S.push_back(  0.094514, 0, 14.330262 ); 
S.push_back(         0, 0, 9.315692 ); 
S.push_back(  0.005532, 0, 69.696126 ); 
S.push_back(         0, 0, 8.578818 ); 
S.push_back(  0.011038, 0, 35.174598 ); 
S.push_back(         0, 0, 4.061354 ); 
S.push_back(  0.005167, 0, 25.657338 ); 
S.push_back(         0, 0, 4.407387 ); 
S.push_back( -0.093196, 0, 9.584562 ); 
S.push_back(         0, 0, 10.380151 ); 
S.push_back( -0.143758, 0, 10.951245 );
S.push_back(         0, 0, 68.584559 ); 
S.push_back( -0.083657, 0, 23.339162 ); 
S.push_back(         0, 0, 5.609678 ); 
S.push_back(  0.027408, 0, 26.928461 ); 
S.push_back(         0, 0, 23.846696 ); 
S.push_back( -0.003647, 0, 77.431314 );
S.push_back(         0, 0, 24.645739 ); 
S.push_back( -0.002587, 0, 86.529164 ); 
%S.push_back(         0, 0, 65.213599 );
S.push_back_G1( x0, y0, theta0); % close curve
S.plot();
axis equal