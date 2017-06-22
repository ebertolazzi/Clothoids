%=============================================================================%
% evalClothoid:  Compute clothoid parameters along a Clothoid curve           %
%                                                                             %
% USAGE: X,Y,TH,K = evalClothoid( x0, y0, theta0, kappa, dkappa, s ) ;        %
%                                                                             %
% On input:                                                                   %
%                                                                             %
%      x0, y0  = coodinate of initial point                                   %
%      theta0  = orientation (angle) of the clothoid at initial point         %
%      kappa   = curvature at initial point                                   %
%      dkappa  = derivative of curvature respect to arclength                 %
%      s       = vector of curvilinear coordinate where to compute clothoid   %
%                                                                             %
% On output:                                                                  %
%                                                                             %
%      X      = X coordinate of the points of the clothoid at s coordinate    %
%      Y      = Y coordinate of the points of the clothoid at s coordinate    %
%      TH     = angle of the clothoid at s coordinate                         %
%      K      = curvature of the clothoid at s coordinate                     %
%                                                                             %
%=============================================================================%
%                                                                             %
%  Autors: Enrico Bertolazzi and Marco Frego                                  %
%          Department of Industrial Engineering                               %
%          University of Trento                                               %
%          enrico.bertolazzi@unitn.it                                         %
%          m.fregox@gmail.com                                                 %
%                                                                             %
%=============================================================================%
function varargout = evalClothoid( x0, y0, theta0, kappa, dkappa, s )
  npts = length(s) ; % how many points to compute
  if nargout < 2
    error( 'evalClothoid:', 'expected at least 2 output arguments' ) ;
  end
  X = zeros(npts,1) ;
  Y = zeros(npts,1) ;
  for k=1:npts    
    [C,S] = GeneralizedFresnelCS( 1, dkappa*s(k)^2, kappa*s(k), theta0 ) ;
    X = x0 + s(k)*C ;
    Y = y0 + s(k)*S ;
  end
  varargout{1} = X ;
  varargout{2} = Y ;
  if nargout > 2
    varargout{3} = theta0 + s*(kappa+s*(dkappa/2)) ; 
  end
  if nargout > 3
    varargout{4} = kappa+s*dkappa ; 
  end
  if nargout > 4
    error( 'evalClothoid:', 'expected at most 4 output arguments' ) ;
  end
end
