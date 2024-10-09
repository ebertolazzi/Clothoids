%=========================================================================%
%                                                                         %
%  Autors: Enrico Bertolazzi                                              %
%          Department of Industrial Engineering                           %
%          University of Trento                                           %
%          enrico.bertolazzi@unitn.it                                     %
%          m.fregox@gmail.com                                             %
%                                                                         %
%=========================================================================%
% Driver test program to check Clothoids lib                              %
%=========================================================================%

close all;

% check constructors
x     = 0:1:19;
y     = [zeros(1,10),ones(1,10)];
%theta = [0,pi/2,-pi/2,0];
if false
y(9)  = y(9)+0.05;
y(10) = y(10)+0.3;
y(11) = y(11)-0.3;
y(12) = y(12)-0.05;
end
if true
y(8)  = y(8)+0.01;
y(9)  = y(9)+0.09;
y(10) = y(10)+0.32;
y(11) = y(11)-0.32;
y(12) = y(12)-0.09;
y(13) = y(13)-0.01;
end

figure('Position',[ 1 1 800 800]);

S = ClothoidList();
S.build_G1( x, y )%, theta );

subplot(5,2,1);
S.plot(400,{'Color','blue','LineWidth',3},{'Color','red','LineWidth',3});
axis equal;

subplot(5,2,2);

fmt1 = {'Color','blue','Linewidth',2};
fmt2 = {'Color','red','Linewidth',2};
S.plotCurvature( 400, fmt1, fmt2 );

for k=[3,5,7,9]

  [ok,maxdK] = S.smooth_quasi_G2(10,1e-10);

  subplot(5,2,k);
  S.plot(400,{'Color','blue','LineWidth',3},{'Color','red','LineWidth',3});
  axis equal;

  subplot(5,2,k+1);
  S.plotCurvature( 400, fmt1, fmt2 );
end


