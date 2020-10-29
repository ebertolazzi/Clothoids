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

SC = ClothoidSplineG2();
%S.ipopt(true);

subplot(3,1,1);

S = SC.buildP2( X, Y );
S.rotate(pi/4,0,0);
len = S.length();

fmt1 = {'Color','red','LineWidth',4};
fmt2 = {'Color','blue','LineWidth',4};
S.plot( 100, fmt1, fmt2);

S1 = ClothoidList();
S1.copy(S);

fmt3 = {'Color','black','LineWidth',3};
fmt4 = {'Color','green','LineWidth',3};
s0 = 0.1*len;
s1 = 0.9*len;
S1.trim( s0, s1 );
S1.plot( 100, fmt3, fmt4);

axis equal

subplot(3,1,2);
S.plotAngle( 100, fmt1, fmt2);
hold on;
S1.plotAngle( 100, fmt3, fmt4);

subplot(3,1,3);
S.plotCurvature( 100, fmt1, fmt2);
hold on;
S1.plotCurvature( 100, fmt3, fmt4);

ds = 0.234;
ss = [0,5];

[x1,y1,theta1,kappa1] = S.evaluate( s0+ds+ss );
[x2,y2,theta2,kappa2] = S1.evaluate( 0+ds+ss );

length(x1)
length(x2)

disp('-----------');
for k=1:length(ss)
  fprintf('x = %g y = %g theta = %g kappa = %g\n', x1(k), y1(k), theta1(k), kappa1(k) );
  fprintf('x = %g y = %g theta = %g kappa = %g\n', x2(k), y2(k), theta2(k), kappa2(k) );
  fprintf('dx = %g dy = %g dtheta = %g dkappa = %g\n', x1(k)-x2(k), y1(k)-y2(k), theta1(k)-theta2(k), kappa1(k)-kappa2(k) );
end

disp('-----------');

ss = [0,1,2,3,4,5];
for k=1:length(ss)
  [x1,y1,theta1,kappa1] = S.evaluate( s0+ds+ss(k) );
  [x2,y2,theta2,kappa2] = S1.evaluate( 0+ds+ss(k) );
  fprintf('x = %g y = %g theta = %g kappa = %g\n', x1, y1, theta1, kappa1 );
  fprintf('x = %g y = %g theta = %g kappa = %g\n', x2, y2, theta2, kappa2 );
  fprintf('dx = %g dy = %g dtheta = %g dkappa = %g\n', x1-x2, y1-y2, theta1-theta2, kappa1-kappa2 );
end