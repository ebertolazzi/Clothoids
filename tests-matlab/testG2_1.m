addpath('../matlab');

%=============================================================================%
%                                                                             %
%  Autors: Enrico Bertolazzi                                                  %
%          Department of Industrial Engineering                               %
%          University of Trento                                               %
%          enrico.bertolazzi@unitn.it                                         %
%          m.fregox@gmail.com                                                 %
%                                                                             %
%=============================================================================%
% Driver test program to check bounding box on clothoid                       %
%=============================================================================%

close all ;
format long ;
x0 = -1 ;
y0 = 0 ;
x1 = 1 ;
y1 = 0 ;

theta0 = pi/3;
kappa0 = [0,1,2,5,10];
theta1 = pi/10;
kappa1 = [0,5,10,50,200];

close all ;
figure(1);

for k=1:4
  %subaxis(2,2,k, 'Spacing', 0.02, 'Padding', 0.04, 'Margin', 0.02);
  %axis tight ;
  subplot(2,2,k);
  for dmax=[0.3,0.6,1]
      dmax
    [ S0, S1, SM, SG, iter ] = buildClothoid3arcG2(x0,y0,theta0,kappa0(k),x1,y1,theta1,kappa1(k),pi,pi/4*dmax) ;
    if iter >= 0
      draw3curve( S0, S1, SM, dmax ~= 1 );
      hold on
      axis equal;
      [X,Y] = pointsOnClothoid( SG, 0:SG.L/100:SG.L );
      plot(X,Y,'-.g','Linewidth',1);
    else
      disp('errore') ;
    end
  end
  title( sprintf('$\\kappa_0=%g,\\quad \\kappa_1=%g%s$',kappa0(k),kappa1(k)),...
         'Interpreter','none'); 
  set(gca,'DataAspectRatio',[1,1,1]) ;
  if false
    set(gca,'XTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]) ;
    set(gca,'XTickLabel',{'$-2$','$-1.5$','$-1$','$-0.5$','$0$','$0.5$','$1$','$1.5$','$2$'}) ;
    set(gca,'YTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]) ;
    set(gca,'YTickLabel',{'$-2$','$-1.5$','$-1$','$-0.5$','$0$','$0.5$','$1$','$1.5$','$2$'}) ;
  end
end

if false
  addpath('tools');
  addpath('tools/matlab2tikz/src');
  matlab2tikz('figure1.tex', ...
              'standalone',true, ...
              'extraaxisoptions',{'title style={font=\Large}','xlabel style={font=\LARGE}','ylabel style={font=\LARGE}','ticklabel style={font=\LARGE}'}, ...
              'extraTikzpictureOptions',{'cap=round'}); 
end

