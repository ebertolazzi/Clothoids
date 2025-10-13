% Author      : Frank E. Curtis
% Description : Driver class for Pipal algorithm.

% Pipal class
classdef Pipal
  
  % Class properties (private access)
  properties (SetAccess = private, GetAccess = private)
    
    i % Input object
    o % Output object
    c % Counter object
    p % Parameter object
    z % Iterate object
    d % Direction object
    a % Acceptance object
    
  end
  
  % Class methods
  methods
    
    % Constructor
    function P = Pipal(name, f_orig, c_orig, g_orig, J_orig, H_orig, x0, bl, bu, l, cl, cu, outfile, algorithm)
      
      % Construct classes
      P.p = Parameter(algorithm);
      P.i = Input(P.p, name, f_orig, c_orig, g_orig, J_orig, H_orig, x0, bl, bu, l, cl, cu);
      P.o = Output(P.i,outfile) ;
      P.c = Counter             ;
      P.z = Iterate(P.p,P.i,P.c);
      P.d = Direction           ;
      P.a = Acceptance          ;
      
    end
    
    % Gets primal solution
    function sol = getSolution(P)
      sol = P.z.getSolution(P.i);
    end
    
    % Optimization algorithm
    function optimize(P)
      
      % Print header and line break
      P.o.printHeader(P.i,    P.z);
      P.o.printBreak (    P.c    );
      
      % Iteration loop
      while ~P.z.checkTermination(P.p,P.i,P.c)
        P.o.printIterate   (        P.c,P.z        );
        P.d.evalStep       (P.p,P.i,P.c,P.z,    P.a);
        P.o.printDirection (            P.z,P.d    );
        P.a.lineSearch     (P.p,P.i,P.c,P.z,P.d    );
        P.o.printAcceptance(                    P.a);
        P.z.updateIterate  (P.p,P.i,P.c,    P.d,P.a);
        P.c.incrementIterationCount                 ;
        P.o.printBreak     (        P.c            );
      end
      
      % Print footer and terminate
      P.o.printFooter(P.p,P.i,P.c,P.z);
      P.o.terminate                   ;
      
    end
    
  end
  
end
