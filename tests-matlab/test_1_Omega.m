%=========================================================================%
%                                                                         %
%  Autors: Enrico Bertolazzi                                              %
%          Department of Industrial Engineering                           %
%          University of Trento                                           %
%          enrico.bertolazzi@unitn.it                                     %
%                                                                         %
%=========================================================================%
% Driver test program to check clothoid computation                       %
%=========================================================================%

addpath('../matlab');
%addpath('../../ipopt/precompiled_mex');

tol = 1e-10;

X = [2.9265642,2.6734362,2.5109322,1.9078122,1.1859282,1.9249962, ...
     2.8265562,0.00468420000000025,-2.826567,-1.9437558,-1.1859438, ...
     -1.9062558,-2.501565,-2.6734386,-2.9265642,-2.6187522,-1.1406318, ...
     -0.8968758,-1.4562558,-1.9062558,-0.00468780000000013,1.9078122, ...
     1.4468682,0.8968722,1.1406282,2.6187522, 2.9265642 ];
Y = [-1.707808758,-1.707808758,-2.367185958,-2.582810358,-2.582810358, ...
     -1.167184758,0.915619242,3.178123242,0.915619242,-1.150000758, ...
     -2.582810358,-2.582810358,-2.393750358,-1.707808758,-1.707808758, ...
     -3.178123242,-3.178123242,-2.989063158,-0.915616758,0.925003242, ...
     2.953123242,0.925003242,-0.915616758,-2.989063158,-3.178123242,-3.178123242, -1.707808758 ];

close all;

S  = ClothoidSplineG2();
S.ipopt(true);
%S.ipopt_check(true);

SPL = cell(8,1);

SPL{1} = S.buildP1( X, Y, -pi, -pi ); %3.141592653589793, 7.647610983392527 );
SPL{3} = S.buildP3( X, Y, SPL{1}.thetaBegin(), SPL{1}.kappaBegin() );

SPL{1} = S.buildP1( X, Y, -pi/2, 0 ); %3.141592653589793, 7.647610983392527 );
SPL{2} = S.buildP2( X, Y );
SPL{4} = S.buildP4( X, Y );
SPL{5} = S.buildP5( X, Y );
SPL{6} = S.buildP6( X, Y );
SPL{7} = S.buildP7( X, Y );
SPL{8} = S.buildP8( X, Y );
SPL{9} = S.buildP9( X, Y );


aa = 0.04;
bb = 1/3-2*aa;

figure('Position',[ 1 1 800 800]);

for k=1:9
  switch(k)
  case 1; subplot('Position',[aa     aa     bb bb]);
  case 2; subplot('Position',[aa+1/3 aa     bb bb]);
  case 3; subplot('Position',[aa+2/3 aa     bb bb]);
  case 4; subplot('Position',[aa     aa+1/3 bb bb]);
  case 5; subplot('Position',[aa+1/3 aa+1/3 bb bb]);
  case 6; subplot('Position',[aa+2/3 aa+1/3 bb bb]);
  case 7; subplot('Position',[aa     aa+2/3 bb bb]);
  case 8; subplot('Position',[aa+1/3 aa+2/3 bb bb]);
  case 9; subplot('Position',[aa+2/3 aa+2/3 bb bb]);
  end
  
  % subplot(3,3,k);
  SPL{k}.plot();
  SPL{k}.plotNormal(0.25,0.25);

  if k == 3
    hold on
    SPL{2}.plot(1000,{'Color','red','LineWidth',1},{'Color','blue','LineWidth',1});
  end
  
  title(sprintf('test P%d',k));
  
  axis equal;
end
