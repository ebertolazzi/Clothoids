function biarc_plot(arc1,arc2,plot_poly)
  breaks1 = fnbrk(arc1,'b');
  fnplt(arc1,[breaks1(1),breaks1(end)],'b',1);
  breaks2 = fnbrk(arc2,'b');
  fnplt(arc2,[breaks2(1),breaks2(end)],'r',1);
  if plot_poly
    %vd = fntlr(bs,2,breaks1);
    %quiver(vd(1,:),vd(2,:),vd(4,:),-vd(3,:));
    xx = arc1.coefs(1,:)./arc1.coefs(3,:);
    yy = arc1.coefs(2,:)./arc1.coefs(3,:);
    plot(xx,yy,':ok');
    xx = arc2.coefs(1,:)./arc2.coefs(3,:);
    yy = arc2.coefs(2,:)./arc2.coefs(3,:);
    plot(xx,yy,':ok');
  end
end