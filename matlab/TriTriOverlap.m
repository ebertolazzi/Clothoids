%======================================================================%
%  TriTriOverlap:  Check if two triangles overlap                      %
%                                                                      %
%  USAGE:                                                              %
%    intersect = TriTriOverlap( p0, p1, p2, q0, q1, q2 );             %
%                                                                      %
%  On input:                                                           %
%                                                                      %
%    p0, p1, p2 = coodinates of first triangle                         %
%    q0, q1, q2 = coodinates of second triangle                        %
%                                                                      %
%  On output:                                                          %
%                                                                      %
%    intersect = true if the triangle overlap                          %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function intersect = TriTriOverlap( p0, p1, p2, q0, q1, q2 )
  error('this function is mex only. Run CompileLib.m script to build')
end