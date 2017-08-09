x0 = -1;
y0 = 0;
x1 = 1;
y1 = 0;
kk = 1;

thmin = 0;
thmax = 0.9*pi;

for theta1=thmin:(thmax-thmin)/3:thmax

  subplot(2,2,kk);

  hold off
  plot([x0,x1],[y0,y1],'k');
  xlim([-1.5,1.5]);
  ylim([-1.5,1.5]);
  axis([-1.5,1.5,-1.5,1.5]);
  hold on

  for theta0=[-pi:pi/20:pi]
    [arc1,arc2,ok] = biarc(x0,y0,theta0,x1,y1,theta1);
    if ok
      biarc_plot(arc1,arc2,false);
    end
  end
  axis equal

  kk = kk+1;
end