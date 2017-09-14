addpath('tools');
addpath('tools/matlab2tikz/src');

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
  subaxis(2,2,k, 'Spacing', 0.02, 'Padding', 0.04, 'Margin', 0.02);
  axis tight ;
  %subplot(2,2,k)
  for thmax0=pi*[0.25,1,1.75]
    [ S0, S1, SM, SG, iter ] = buildClothoid3arcG2(x0,y0,theta0,kappa0(k),x1,y1,theta1,kappa1(k),thmax0,thmax0) ;
    if iter >= 0
      draw3curve( S0, S1, SM, thmax0 ~= pi );
      hold on
      axis equal;
      title( sprintf('k0=%g, k1=%g',kappa0(k),kappa1(k))); 
    else
      disp('errore') ;
    end
  end
  set(gca,'DataAspectRatio',[1,1,1]) ;
  if k == 4
    set(gca,'XTick',[-2,-1,0,1,2]) ;
    set(gca,'XTickLabel',{'$-2$','$-1$','$0$','$1$','$2$'}) ;
    set(gca,'YTick',[-2,-1,0,1,2,3]) ;
    set(gca,'YTickLabel',{'$-2$','$-1$','$0$','$1$','$2$','$3$'}) ;
  else
    set(gca,'XTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]) ;
    set(gca,'XTickLabel',{'$-2$','$-1.5$','$-1$','$-0.5$','$0$','$0.5$','$1$','$1.5$','$2$'}) ;
    set(gca,'YTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]) ;
    set(gca,'YTickLabel',{'$-2$','$-1.5$','$-1$','$-0.5$','$0$','$0.5$','$1$','$1.5$','$2$'}) ;
  end
end

if true
  matlab2tikz('figure1.tex', ...
              'standalone',true, ...
              'extraaxisoptions',{'xlabel style={font=\LARGE}','ylabel style={font=\LARGE}','ticklabel style={font=\LARGE}'}, ...
              'extraTikzpictureOptions',{'cap=round'}); 
end

