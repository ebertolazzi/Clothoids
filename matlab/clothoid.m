%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef clothoid < handle
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    %% Constructor - Create a new C++ class instance 
    function this = clothoid( x0, y0, theta0, x1, y1, theta1 )
      this.objectHandle = ClothoidMexWrapper('new', x0, y0, theta0, x1, y1, theta1 );
    end
        
    %% Destructor - Destroy the C++ class instance
    function delete(this)
      ClothoidMexWrapper('delete', this.objectHandle );
    end

    % XY, [X,Y]
    function varargout = eval(this, s)
      varargout = ClothoidMexWrapper('eval', this.objectHandle, s );
    end

    function [x0,y0,theta0,kappa0,dk,L] = getPars(this)
      [x0,y0,theta0,kappa0,dk] = ClothoidMexWrapper('getPars', this.objectHandle );
    end
    
    % add......

    function info(this)
      ClothoidMexWrapper('info', this.objectHandle);
    end
  
  end
end