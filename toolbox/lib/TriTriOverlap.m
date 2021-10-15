%> TriTriOverlap:  Check if two triangles overlap
%>
%> **Usage:**
%>
%>     intersect = TriTriOverlap( p0, p1, p2, q0, q1, q2 );
%>
%> **On input:**
%>
%> - p0, p1, p2 = coordinates of first triangle
%> - q0, q1, q2 = coordinates of second triangle
%>
%> **On output:**
%>
%> - intersect = true if the triangle overlap
%>
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function intersect = TriTriOverlap( p0, p1, p2, q0, q1, q2 )
  error('this function is mex only. Run CompileClothoidsLib.m script to build')
end