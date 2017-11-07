function [arc1,arc2,ok] = biarc_minlen(x0,y0,th0,x1,y1,th1)
  [arc1,arc2,ok] = biarc(x0,y0,th0,x1,y1,th1);
  %mlen = arc1.length^2+arc2.length^2;
  mlen = abs(arc1.curvature)+abs(arc2.curvature);
  for th=-pi:pi/1000:pi
    [a1,a2,ok] = biarc(x0,y0,th0,x1,y1,th1,th);
    if ok
      %len = a1.length^2+a2.length^2;
      len = abs(a1.curvature)+abs(a2.curvature);
      if len < mlen
        arc1 = a1 ;
        arc2 = a2 ;
      end
    end
  end
end