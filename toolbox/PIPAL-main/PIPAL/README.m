% README.m
%
% Please cite:
%   Frank E. Curtis.  "A Penalty-Interior-Point Algorithm for Nonlinear Constrained
%   Optimization."  Mathematical Programming Computation, 4(2):181–209, 2012.
%
% This code offers penalty-interior-point algorithms for solving nonlinear
% constrained optimization problems of the form
%   minimize f(x) subject to cE(x) = 0, cI(x) <= 0,
% where f, cE, and cI are assumed to be continuously differentiable in R^n.
% The code solves problems in the form of AMPL .nl files *ONLY*.
%
% Solving an AMPL problem instance requires the following input.
%   name      ~ (string) name of problem
%   nl        ~ (string) name of AMPL .nl file
%   outfile   ~ (string) name of output file
%   algorithm ~ (integer in {0,1}) algorithm number
%               (0 ~ conservative algorithm, 1 ~ adaptive algorithm)
% The problem is solved and the solution obtained via the commands
%   >> P = Pipal(name,nl,outfile,algorithm);
%   >> P.optimize;
%   >> s = P.getSolution;
% The solution structure s yields the following output.
%   s.x ~ final primal point
%   s.l ~ final dual point
%
% (An amplfunc.mex file that reads AMPL .nl files is necessary; please see
% D. Gay, ``Hooking Your Solver to AMPL,'' Technical Report, Computing
% Sciences Research Center, Bell Laboratories, Murray Hill, NJ, USA.)