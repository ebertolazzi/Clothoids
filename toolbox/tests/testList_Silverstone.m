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

S = ClothoidList();
S.push_back( 0, 0, 0.079964, 0, 0, 212.092029 ); 
S.push_back( -0.012331, 0, 109.618459 ); 
S.push_back(         0, 0, 121.126141 ); 
S.push_back( -0.001926, 0, 129.688509 ); 
S.push_back(         0, 0, 176.642769 ); 
S.push_back(  0.005072, 0, 104.112724 ); 
S.push_back(         0, 0, 6.544815 ); 
S.push_back( -0.025013, 0, 37.020151 ); 
S.push_back(         0, 0, 61.700509 ); 
S.push_back(  0.005864, 0, 85.660725 ); 
S.push_back(         0, 0, 5.232855 ); 
S.push_back(  0.014306, 0, 35.558911 ); 
S.push_back(         0, 0, 37.161532 ); 
S.push_back( -0.011938, 0, 64.163520 ); 
S.push_back(         0, 0, 3.459784 ); 
S.push_back( -0.009904, 0, 95.218635 ); 
S.push_back(         0, 0, 32.800804 ); 
S.push_back(  0.006176, 0, 89.902292 ); 
S.push_back(         0, 0, 578.623657 ); 
S.push_back( -0.001759, 0, 156.859397 ); 
S.push_back(         0, 0, 11.870129 ); 
S.push_back( -0.011746, 0, 83.076858 ); 
S.push_back(         0, 0, 12.699057 ); 
S.push_back( -0.013701, 0, 69.926391 ); 
S.push_back(         0, 0, 29.777322 ); 
S.push_back(  0.002070, 0, 140.586962 ); 
S.push_back(         0, 0, 184.860552 ); 
S.push_back(  0.054896, 0, 11.582006 );
S.push_back(         0, 0, 1.819584 ); 
S.push_back(  0.047439, 0, 20.573067 ); 
S.push_back(         0, 0, 28.752633 ); 
S.push_back( -0.024724, 0, 33.441840 ); 
S.push_back(         0, 0, 3.963619 ); 
S.push_back( -0.015993, 0, 44.528072 ); 
S.push_back(         0, 0, 7.951332 ); 
S.push_back( -0.009115, 0, 83.819087 ); 
S.push_back(         0, 0, 5.849496 ); 
S.push_back( -0.007803, 0, 88.137874 ); 
S.push_back(         0, 0, 491.418500 ); 
S.push_back(  0.034222, 0, 46.113541 ); 
S.push_back(         0, 0, 46.264337 ); 
S.push_back( -0.031383, 0, 27.345066 ); 
S.push_back(         0, 0, 203.645328 ); 
S.push_back( -0.012461, 0, 85.722960 );
S.push_back(         0, 0, 139.761254 ); 
S.push_back(  0.016207, 0, 55.369995 ); 
S.push_back(         0, 0, 4.479479 ); 
S.push_back(  0.013413, 0, 53.028182 ); 
S.push_back(         0, 0, 105.388975 ); 
S.push_back(  0.027805, 0, 43.375087 ); 
S.push_back(         0, 0, 5.499410 ); 
S.push_back(  0.033446, 0, 36.462857 ); 
S.push_back(         0, 0, 74.222954 ); 
S.push_back( -0.020038, 0, 89.813675 ); 
S.push_back(         0, 0, 8.514349 ); 
S.push_back( -0.022536 , 0, 74.953886 );
S.push_back(         0, 0, 126.296448 ); 
S.push_back( -0.004020, 0, 112.003093 ); 
S.push_back(         0, 0, 13.793995 ); 
S.push_back( -0.006581, 0, 86.067741 ); 
%S.push_back(         0, 0, 207.629231 ); 
S.push_back_G1( 0, 0, 0.079964 ); % close curve
S.plot();
S.save('silverstone.txt');
axis equal
