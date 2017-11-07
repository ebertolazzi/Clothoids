classdef ClothoidCurve < handle
    % clothoid: MATLAB class wrapper to the underlying C++ class
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    
    methods
        %% Constructor - Create a new C++ class instance
        function this = ClothoidCurve( varargin )
            % clothoid: constructor for the clothoid object.
            % Usage:
            %    ref = clothoid()
            %    ref = clothoid( x0, y0, theta0, k0, dk, L )
            %    ref = clothoid( x0, y0, theta0, k0, dk, smin, smax )
            %
            % On input:
            %    x0, y0: coordinate of initial point
            %    theta0: orientation of the clothoid at initial point
            %    k0:     curvature of the clothoid at initial point
            %    dk:     derivative of curvature respect to arclength
            %    L:      length of curve from initial to final point
            %    smin:   initial curvilinear coordinate of the curve
            %    smax:   final curvilinear coordinate of the curve
            %
            %  On output:
            %    ref: reference handle to the object instance
            
            this.objectHandle = ClothoidMexWrapper('new', varargin{:} );
        end
        
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            ClothoidMexWrapper('delete', this.objectHandle );
        end
        
        
        %% Build
        function build( this, varargin )
            % build: method to build the clothoid from known parameters
            % Usage:
            %    ref.build( x0, y0, theta0, k0, dk, L )
            %    ref.build( x0, y0, theta0, k0, dk, smin, smax )
            %
            % On input:
            %    x0, y0: coordinate of initial point
            %    theta0: orientation of the clothoid at initial point
            %    k0:     curvature of the clothoid at initial point
            %    dk:     derivative of curvature respect to arclength
            %    L :     length of curve from initial to final point
            %    smin:   initial curvilinear coordinate of the curve
            %    smax:   final curvilinear coordinate of the curve
            
            ClothoidMexWrapper('build', this.objectHandle, varargin{:} );
        end
        
        function build_G1( this, x0, y0, theta0, x1, y1, theta1 )
            % build_G1: method to build the interpolating G1 clothoid arc
            % Usage:
            %    ref.build_G1( x0, y0, theta0, x1, y1, theta1 )
            %
            % On input:
            %    x0, y0: coordinate of initial point
            %    theta0: orientation of the clothoid at initial point
            %    x1, y1: coordinate of final point
            %    theta1: orientation of the clothoid at final point
            
            ClothoidMexWrapper('build_G1', this.objectHandle, ...
                x0, y0, theta0, x1, y1, theta1 );
        end
        
        function res = build_forward( this, x0, y0, theta0, k0, x1, y1 )
            % build_forward: method to build the interpolating clothoid arc
            % Usage:
            %    res = ref.build_forward( x0, y0, theta0, k0, x1, y1 )
            %
            % On input:
            %    x0, y0: coordinate of initial point
            %    theta0: orientation of the clothoid at initial point
            %    k0:     curvature of the clothoid at initial point
            %    x1, y1: coordinate of final point
            %
            % On output:
            %    res: true iff the interpolation was successful
            
            res = ...
                ClothoidMexWrapper('build_forward', this.objectHandle, ...
                                   x0, y0, theta0, k0, x1, y1 );
        end
        
        
        %% Eval
        function varargout = eval(this, s)
            % eval: method to eval the curve at curvilinear abscissa s
            % Usage:
            %    [x,y] = ref.eval( s )
            %    [x,y,theta,kappa] = ref.eval( s )
            %    
            % On input:
            %    s: curvilinear coordinates where to evaluate the curve
            %       (scalar or vector)
            %
            % On output:
            %    x, y:  coordinates of the curve 
            %    theta: orientation of the curve
            %    kappa: curvature of the curve
            if nargout == 1
              [x,y] = ClothoidMexWrapper('eval', this.objectHandle, s );
              varargout{1} = [x,y] ;
            else
              [x,y,varargout{2},varargout{3}] = ClothoidMexWrapper('eval', this.objectHandle, s );
              varargout{1} = [x,y] ;
            end
            
        end

        
        %%
        function [X,Y,S] = points(this, ds)
          [~,~,~,~,~,smin,smax] = ClothoidMexWrapper('getPars', this.objectHandle );
          S = smin:(smax-smin)/ds:smax ;
          [X,Y] = ClothoidMexWrapper('eval', this.objectHandle, S );
        end
 
        function [x0,y0,theta0,k0,dk,smin,smax] = getPars(this)
            % getPars: method to get the parameters of the curve
            % Usage:
            %    [x0,y0,theta0,k0,dk,smin,smax] = ref.getPars()
            %    
            % On output:
            %    x0, y0: coordinate of initial point
            %    theta0: orientation of the clothoid at initial point
            %    k0:     curvature of the clothoid at initial point
            %    dk:     derivative of curvature respect to arclength
            %    smin:   initial curvilinear coordinate of the curve
            %    smax:   final curvilinear coordinate of the curve
            
            [x0,y0,theta0,k0,dk,smin,smax] = ...
                ClothoidMexWrapper('getPars', this.objectHandle );
        end
        
 
        function len = length(this)
            [~,~,~,~,~,smin,smax] = ...
                ClothoidMexWrapper('getPars', this.objectHandle );
            len = smax-smin ;
        end
        
        %% Update
        function trim(this, smin, smax)
            % trim: method to trim the clothoid curve
            % Usage:
            %    ref.trim(smin, smax)
            %    
            % On input:
            %    smin:   initial curvilinear coordinate of the curve
            %    smax:   final curvilinear coordinate of the curve
            
            ClothoidMexWrapper('trim', this.objectHandle, smin, smax );
        end
        
        function changeOrigin(this, s0)
            % changeOrigin: method to change the origin of the clothoid 
            %               curve to s0
            % Usage:
            %    ref.changeOrigin(s0)
            %    
            % On input:
            %    s0: curvilinear coordinate of the origin of the new curve
            
            ClothoidMexWrapper('changeOrigin', this.objectHandle, s0 );
        end
        
        function rotate(this, angle, cx, cy)
            % rotate: method to rotate the clothoid curve
            % Usage:
            %    ref.rotate(angle, cx, cy)
            %    
            % On input:
            %    angle: the angle of rotation
            %    cx, cy: coordinates of the centre of rotation
            
            ClothoidMexWrapper('rotate', this.objectHandle, ...
                               angle, cx, cy );
        end
        
        function translate(this, tx, ty)
            % translate: method to translate the clothoid curve
            % Usage:
            %    ref.translate(tx, ty)
            %    
            % On input:
            %    tx, ty: horizontal and vertical translation
            
            ClothoidMexWrapper('translate', this.objectHandle, tx, ty );
        end

        function moveOrigin(this, newX0, newY0)
            % moveOrigin: method to move the origin of the clothoid curve 
            % Usage:
            %    ref.moveOrigin(newX0, newY0)
            %    
            % On input:
            %    newX0, newY0: new coordinates of initial point
            
            ClothoidMexWrapper('moveOrigin', this.objectHandle, ...
                               newX0, newY0 );
        end
        
        function scale(this, s)
            % scale: method to scale the clothoid curve 
            % Usage:
            %    ref.scale(newX0, newY0)
            %    
            % On input:
            %    newX0, newY0: new coordinates of initial point
            
            ClothoidMexWrapper('scale', this.objectHandle, s );
        end
        
        function reverse(this)
            % reverse: method to reverse the clothoid curve 
            % Usage:
            %    ref.reverse()
            
            ClothoidMexWrapper('reverse', this.objectHandle );
        end
        
        
        %% Utils
        function lineH = plot(this, step, varargin)
            % plot: method to plot the clothoid curve
            % Usage:
            %    lineH = ref.plot()
            %    lineH = ref.plot(0.2)
            %    lineH = ref.plot(0.2, 'r', 'LineWidth', 3)
            %    
            % On input:
            %    step:     the distance between consecutive samples used to
            %              draw the curve
            %    varargin: optional arguments passed to the MATLAB plot
            %              function
            %
            % On output:
            %    lineH: the handle to the line object 
            %           (i.e. the output of the MATLAB plot function)
            
            if nargin<2
                step = 0.1;
            end
            [x0, y0, theta0, k0, dk, smin, smax] = this.getPars();
            XY = pointsOnClothoid(x0, y0, theta0, k0, dk, ...
                [smin:step:smax smax]);
            lineH = plot(XY(1,:),XY(2,:), varargin{:});
        end
        
        
    end
end