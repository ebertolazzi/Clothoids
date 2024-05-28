%-------------------------------------------------------+
%                                                       |
% Copyright (C) 2022                                    |
%                                                       |
%        , __                 , __                      |
%       /|/  \               /|/  \                     |
%        | __/ _   ,_         | __/ _   ,_              |
%        |   \|/  /  |  |   | |   \|/  /  |  |   |      |
%        |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/     |
%                          /|                   /|      |
%                          \|                   \|      |
%                                                       |
%     Enrico Bertolazzi                                 |
%     Dipartimento di Ingegneria Industriale            |
%     Universita` degli Studi di Trento                 |
%     email: enrico.bertolazzi@unitn.it                 |
%                                                       |
%-------------------------------------------------------+
%
%      _    _           _____ _  _    ___
%     / \  | | __ _  __|___  | || |  ( _ )
%    / _ \ | |/ _` |/ _ \ / /| || |_ / _ \
%   / ___ \| | (_| | (_) / / |__   _| (_) |
%  /_/   \_\_|\__, |\___/_/     |_|  \___/
%             |___/
%
%
% Based on the paper:
%
% - **G. E. Alefeld, Florian A Potra, Yixun Shi**,
%   *Algorithm 748: enclosing zeros of continuous functions*,
%   ACM Transactions on Mathematical Software, vol 21, N.3, 1995
%
classdef Algo748 < matlab.mixin.Copyable

  properties (SetAccess = private, Hidden = true)
    m_mu;
    m_fun;
    m_num_fun_eval;
    m_num_iter_done;
    m_tolerance;
  end

  methods (Access = private, Hidden = true)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function fx = f_evaluate( self, x )
      self.m_num_fun_eval = self.m_num_fun_eval+1;
      fx = feval( self.m_fun, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    [a,b,c,d,fa,fb,fc,fd] = bracketing( self, P, FP )
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    res = pzero( ~, P, FP )
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    res = newton_quadratic( ~, niter, P, FP )
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function set_tolerance( self, B )
      self.m_tolerance = 2*(eps + 2*abs(B)*eps);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function yes_no = all_different( ~, V )
      yes_no = false;
      for i=1:length(V)-1
        for j=i+1:length(V)
          if V(i) == V(j)
            return;
          end
        end
      end
      yes_no = true;
    end
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Algo748()
      self.m_num_iter_done = 0;
      self.m_num_fun_eval  = 0;
      self.m_mu            = 0.5;
      self.m_tolerance     = 1e-8;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = num_fun_eval( self )
      res = self.m_num_fun_eval;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = num_iter_done( self )
      res = self.m_num_iter_done;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    [sol,iter] = eval( self, a, b, fun )
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
