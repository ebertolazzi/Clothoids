addpath('tools');
addpath('tools/matlab2tikz/src');

close all ;
format long ;
x0 = -1 ;
y0 = 0 ;
x1 = 1 ;
y1 = 0 ;

theta0 = pi/3;
kappa0 = 0;
theta1 = [-pi/2,pi/2,-pi/2,pi/2,-pi/2,pi/2,];% [-pi/2,-pi/8,0,pi/4,pi/2];
kappa1 = [0,0,-4,-4,4,4];

close all ;
figure('Position',[1,1,900,800]);

AXEL = [-1 1 -0.25 1.25] ;
AXER = [-1 1 -0.75 0.75] ;

for k=1:6
  subaxis(3,2,k, 'Spacing', 0.01, 'Padding', 0.02, 'Margin', 0.02);
  axis tight ;
  
  [ S0, S1, SM, SG, iter ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,x1,y1,theta1(k),kappa1(k)) ;
  if iter >= 0
    draw3curve( S0, S1, SM, false );
    hold on
    if mod(k,2) == 0
      axis(AXER);
    else
      axis(AXEL);
    end
    axis equal;
    title( sprintf('kappa1 = %g',kappa1(k))) ; 

    [X,Y] = pointsOnClothoid( SG, 0:SG.L/100:SG.L );
     plot(X,Y,'-m','Linewidth',1);

    %[arc1,arc2,ok] = biarc(x0,y0,theta0,x1,y1,theta1(k));
    %if ok
    %  biarc_plot(arc1,arc2,false);
    %end
  else
    disp('errore') ;
  end
  set(gca,'DataAspectRatio',[1,1,1]) ;
  set(gca,'XTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]) ;
  set(gca,'XTickLabel',{'$-2$','$-1.5$','$-1$','$-0.5$','$0$','$0.5$','$1$','$1.5$','$2$'}) ;
  set(gca,'YTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]) ;
  set(gca,'YTickLabel',{'$-2$','$-1.5$','$-1$','$-0.5$','$0$','$0.5$','$1$','$1.5$','$2$'}) ;
end

if true
  matlab2tikz('figure2.tex', ...
              'standalone',true, ...
              'extraaxisoptions',{'title style={font=\Large}','xlabel style={font=\LARGE}','ylabel style={font=\LARGE}','ticklabel style={font=\LARGE}'}, ...
              'extraTikzpictureOptions',{'cap=round'}); 
end
