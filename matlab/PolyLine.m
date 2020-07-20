classdef PolyLine < CurveBase
  %% MATLAB class wrapper for the underlying C++ class
  methods
    function self = PolyLine( )
      self@CurveBase( 'PolyLineMexWrapper' );
      self.objectHandle = PolyLineMexWrapper( 'new' );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function str = is_type( ~ )
      str = 'PolyLine';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function build( self, x, y )
      PolyLineMexWrapper( 'build', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function approx( self, obj, tol )
      PolyLineMexWrapper( 'approx', self.objectHandle, ...
                                    obj.objectHandle,  ...
                                    tol,               ...
                                    obj.is_type() );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x, y ] = polygon( self )
      [ x, y ] = PolyLineMexWrapper( 'polygon', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, varargin )
      if nargin > 1
        fmt1 = varargin{1};
      else
        fmt1 = {'Color','red','LineWidth',3};
      end
      if nargin > 2
        fmt2 = varargin{2};
      else
        fmt2 = {'Color','blue','LineWidth',3};
      end
      [ x, y ] = self.polygon();
      for k=2:length(x)
        if mod(k,2) == 0
          plot( x(k-1:k), y(k-1:k), fmt1{:} );
        else
          plot( x(k-1:k), y(k-1:k), fmt2{:} );
        end
        hold on;
      end
    end
  end
end
