% Author      : Frank E. Curtis
% Description : Class for handling algorithm output.

% Output class
classdef Output < handle
  
  % Class properties (private access)
  properties (SetAccess = private, GetAccess = private)
    
    s % Output stream
    l % Line break string
    q % Quantities string
    n % Last iterate string
    
  end
  
  % Class methods
  methods
    
    % Constructor
    function o = Output(i,outfile)
      
      % Start clock
      tic;
      
      % Set output stream
      o.s = fopen(outfile,'w');
      
      % Assert output stream has been opened
      assert(o.s~=-1,sprintf('PIPAL: Failed to open %s.out.',i.id));
      
      % Store output strings
      o.l = '======+=========================+====================================+=========================+===========================================================================+=======================';
      o.q = 'Iter. |  Objective     Infeas.  |  Pen. Par.   I.P. Par.  Opt. Error |    Merit     P.I.P. Err.|    Shift    ||P.Step||  ||D.Step||   Lin. Red.    Quad. Red.    Quality   | Pri. Step.  Dual Step.';
      o.n =                                                                        '-----------  ---------- | ----------  ----------  ----------  -----------  -----------  ----------- | ----------  ----------';
      
    end
    
    % Print acceptance information
    function printAcceptance(o,a)
      
      % Print steplengths
      fprintf(o.s,'%.4e  %.4e',a.p,a.d);
      if a.s == 1, fprintf(o.s,' SOC'); end;
      fprintf(o.s,'\n');
      
    end
    
    % Print break in output
    function printBreak(o,c)
      
      % Print break every 20 iterations
      if mod(c.k,20) == 0, fprintf(o.s,'%s\n%s\n%s\n',o.l,o.q,o.l); end;
      
    end
    
    % Print direction information
    function printDirection(o,z,d)
      
      % Print direction information
      fprintf(o.s,'%+.4e  %.4e | %.4e  %.4e  %.4e  %+.4e  %+.4e  %+.4e | ',z.phi,z.kkt(3),z.shift,d.x_norm,d.l_norm,d.ltred,d.qtred,d.m);
      
    end
    
    % Print output footer
    function printFooter(o,p,i,c,z)
      
      % Print final iterate information
      o.printIterate(c,z);
      
      % Print close of algorithm output
      fprintf(o.s,'%s\n%s\n\n',o.n,o.l);
      
      % Get solver result
      b = z.checkTermination(p,i,c);
      
      % Print solver result
      fprintf(o.s,'Final result\n');
      fprintf(o.s,'============\n');
      if b == 0, fprintf(o.s,'  EXIT: No termination message set\n'       ); end;
      if b == 1, fprintf(o.s,'  EXIT: Optimal solution found\n'           ); end;
      if b == 2, fprintf(o.s,'  EXIT: Infeasible stationary point found\n'); end;
      if b == 3, fprintf(o.s,'  EXIT: Iteration limit reached\n'          ); end;
      if b == 4, fprintf(o.s,'  EXIT: Invalid bounds\n'                   ); end;
      if b == 5, fprintf(o.s,'  EXIT: Function evaluation error\n'        ); end;
      fprintf(o.s,'\n');
      
      % Print iterate quantities
      fprintf(o.s,'Final values\n');
      fprintf(o.s,'============\n');
      fprintf(o.s,'  Objective function........................ : %+e\n',z.fu);
      fprintf(o.s,'  Feasibility violation..................... : %+e\n',z.vu);
      fprintf(o.s,'  Optimality error (feasibility)............ : %+e\n',z.kkt(1));
      fprintf(o.s,'  Optimality error (penalty)................ : %+e\n',z.kkt(2));
      fprintf(o.s,'  Optimality error (penalty-interior-point). : %+e\n',z.kkt(3));
      fprintf(o.s,'  Penalty parameter......................... : %+e\n',z.rho);
      fprintf(o.s,'  Interior-point parameter.................. : %+e\n',z.mu);
      fprintf(o.s,'\n');
      
      % Print counters
      fprintf(o.s,'Final counters\n');
      fprintf(o.s,'==============\n');
      fprintf(o.s,'  Iterations................................ : %d\n',c.k);
      fprintf(o.s,'  Function evaluations...................... : %d\n',c.f);
      fprintf(o.s,'  Gradient evaluations...................... : %d\n',c.g);
      fprintf(o.s,'  Hessian evaluations....................... : %d\n',c.H);
      fprintf(o.s,'  Matrix factorizations..................... : %d\n',c.M);
      fprintf(o.s,'  CPU seconds............................... : %d\n',ceil(toc));
      
    end
    
    % Print output header
    function printHeader(o,i,z)
      
      % Print problem name
      fprintf(o.s,'Problem name\n');
      fprintf(o.s,'============\n');
      fprintf(o.s,'  %s\n',i.id);
      fprintf(o.s,'\n');
      
      % Print problem size
      fprintf(o.s,'Problem size\n');
      fprintf(o.s,'============\n');
      fprintf(o.s,'  Number of variables....................... : %8d\n',i.nV);
      fprintf(o.s,'  Number of equality constraints............ : %8d\n',i.nE);
      fprintf(o.s,'  Number of inequality constraints.......... : %8d\n',i.nI);
      fprintf(o.s,'\n');
      
      % Print problem sparsity
      fprintf(o.s,'Problem sparsity\n');
      fprintf(o.s,'================\n');
      fprintf(o.s,'  Nonzeros in Hessian of Lagrangian......... : %8d\n',z.Hnnz);
      fprintf(o.s,'  Nonzeros in equality constraint Jacobian.. : %8d\n',z.JEnnz);
      fprintf(o.s,'  Nonzeros in inequality constraint Jacobian : %8d\n',z.JInnz);
      fprintf(o.s,'\n');
      
    end
    
    % Print iterate information
    function printIterate(o,c,z)
      
      % Print iterate information
      fprintf(o.s,'%5d | %+.4e  %.4e | %.4e  %.4e  %.4e | ',c.k,z.f,z.v,z.rho,z.mu,z.kkt(2));
      
    end
    
    % Terminate output
    function terminate(o)
      
      % Close nonstandard output stream
      if ~ismember(o.s,[0 1 2]), fclose(o.s); end;
      
    end
    
  end
  
end
