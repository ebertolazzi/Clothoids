addpath('../matlab');

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

%format long;
%addpath(fullfile(pwd,'matlab2tikz','src'));
%addpath('tools');

close all;

p=[0,1;0,0]; u=[0.5,-0.1;-0.25,0.5];
plot(p(1,:),p(2,:),'k');
hold on
biarc = rscvn(p,u); breaks = fnbrk(biarc,'b');
fnplt(biarc,breaks(1:2),'b',3), fnplt(biarc,breaks(2:3),'r',3);
vd = fntlr(biarc,2,breaks);
quiver(vd(1,:),vd(2,:),vd(4,:),-vd(3,:)),
hold off

title('Matlab Manual');

set(gca,'DataAspectRatio',[1,1,1]);
set(gca,'XTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]);
set(gca,'XTickLabel',{'$-2$','$-1.5$','$-1$','$-0.5$','$0$','$0.5$','$1$','$1.5$','$2$'});
set(gca,'YTick',[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]);
set(gca,'YTickLabel',{'$-2$','$-1.5$','$-1$','$-0.5$','$0$','$0.5$','$1$','$1.5$','$2$'});

% to produce tikz figure
if false
  matlab2tikz('figure0.tex', ...
              'standalone',true, ...
              'extraaxisoptions',{'xlabel style={font=\LARGE}','ylabel style={font=\LARGE}','ticklabel style={font=\LARGE}'}, ...
              'extraTikzpictureOptions',{'cap=round'}); 
end


