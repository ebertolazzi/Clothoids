function bspline_plot(bs,plot_poly)
  breaks = fnbrk(bs,'b');
  fnplt(bs,[breaks(1),breaks(end)],'b',1);
  if plot_poly
    vd = fntlr(bs,2,breaks);
    quiver(vd(1,:),vd(2,:),vd(4,:),-vd(3,:));
    xx = bs.coefs(1,:)./bs.coefs(3,:);
    yy = bs.coefs(2,:)./bs.coefs(3,:);
    plot(xx,yy,':ok');
  end
end