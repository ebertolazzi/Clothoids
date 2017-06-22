%======================================================================%
% pointsOnClothoid:  Compute points on a clothoid curve.               %
%                    Used for plotting purpose.                        %
%                                                                      %
% USAGE: XY    = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;  %
% USAGE: [X,Y] = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;  %
% USAGE: XY    = pointsOnClothoid( x0, y0, theta0, k, dk, L ) ;        %
% USAGE: [X,Y] = pointsOnClothoid( x0, y0, theta0, k, dk, L ) ;        %
% USAGE: XY    = pointsOnClothoid( clot, npts ) ;                      %
% USAGE: [X,Y] = pointsOnClothoid( clot, npts ) ;                      %
%                                                                      %
% On input:                                                            %
%                                                                      %
%  x0, y0  = coodinate of initial point                                %
%  theta0  = orientation (angle) of the clothoid at initial point      %
%  k       = curvature at initial point                                %
%  dk      = derivative of curvature respect to arclength              %
%  L       = the lenght of the clothoid curve or a vector of length    %
%            where to compute the clothoid values                      %
%  npts    = number of points along the clothoid                       %
%                                                                      %
% In alternative                                                       %
%  clot    = structure with the field x0, y0, kappa, dkappa, L         %
%                                                                      %
%                                                                      %
% On output: (1 argument)                                              %
%                                                                      %
%  XY = matrix 2 x NPTS whose column are the points of the clothoid    %
%                                                                      %
% On output: (2 argument)                                              %
%                                                                      %
%   X  = matrix 1 x NPTS X coordinate of points of the clothoid        %
%   Y  = matrix 1 x NPTS Y coordinate of points of the clothoid        %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function varargout = pointsOnClothoid( varargin )
  
  if nargin == 2
    if isstruct(varargin{1})
      x0     = varargin{1}.x0 ;
      y0     = varargin{1}.y0 ;
      theta0 = varargin{1}.theta0 ;
      kappa  = varargin{1}.kappa ;
      dkappa = varargin{1}.dkappa ;
      L      = varargin{1}.L ;
      npts   = varargin{2} ;
      tvec   = [0:L/(npts-1):L] ;
    else
      error('expexted struct as first arument') ;
    end
  elseif nargin == 6 || nargin == 7
    x0     = varargin{1} ;
    y0     = varargin{2} ;
    theta0 = varargin{3} ;
    kappa  = varargin{4} ;
    dkappa = varargin{5} ;
    L      = varargin{6} ;
    if nargin == 7
      npts = varargin{7} ;
      tvec = [0:L/(npts-1):L] ;
    else
      tvec = L ;
    end
  else
    error('expexted 2,6 or 7 input arguments') ;
  end
  
  X = [] ;
  Y = [] ;
  for t=tvec
    [C,S] = GeneralizedFresnelCS( 1, dkappa*t^2, kappa*t, theta0 ) ;
    X = [ X x0 + t*C ] ;
    Y = [ Y y0 + t*S ] ;
  end
  if nargout > 1
    varargout{1} = X ;
    varargout{2} = Y ;
  else
    varargout{1} = [X ; Y] ;
  end
end
