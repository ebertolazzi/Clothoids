%======================================================================%
%  intersectClothoid:  Compute intersections betweed clothoids         %
%                                                                      %
%  USAGE:                                                              %
%    s1, s2 = ClothoidCurveIntertsectMexWrapper( clot1, clot2 ) ;      %
%    s1, s2 = ClothoidCurveIntertsectMexWrapper( clot1, clot2,         %
%                                                offs1, offs2 ) ;      %
%                                                                      %
%  On input:                                                           %
%                                                                      %
%    clot1 = object pointer of the first clothoid                      %
%    clot2 = object pointer of the second clothoid                     %
%                                                                      %
%  On output:                                                          %
%                                                                      %
%   s1 = curvilinear coordinates of intersections on clot1             %
%   s2 = curvilinear coordinates of intersections on clot2             %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function [s1, s2] = intersectClothoid( clot1, clot2, varargin ) ;
  [s1, s2] = ClothoidCurveIntertsectMexWrapper( clot1.obj_handle(), ...
                                                clot2.obj_handle(), ...
                                                varargin{:} );
end