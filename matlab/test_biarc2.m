x0=0;
y0=0;
x1=1;
y1=0;

subplot(2,1,1);
hold off
hold on

subplot(2,1,2);
hold off
hold on

theta0=pi/12;
theta1=-pi/4;
a0=theta0+pi/2;
a1=theta1+pi/2;

  p=[x0,x1;y0,y1];
  u=[cos(a0),cos(a1);sin(a0),sin(a1)];

  subplot(2,1,1);
  plot(p(1,:),p(2,:),'k')
  bi_arc = fnrfn(rscvn(p,u),[0.5,1.5]);
  bspline_plot(bi_arc,true);

  subplot(2,1,2);
  plot(p(1,:),p(2,:),'k')
  [arc1,arc2] = biarc(x0,y0,theta0,x1,y1,theta1);
  biarc_plot(arc1,arc2,true);

subplot(2,1,1);
axis equal

subplot(2,1,2);
axis equal
