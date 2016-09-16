classdef ConstPressure < barotropes.Barotrope
    %CONSTPRESSURE A test implementation of constant pressure barotrope.
    
    properties (Access = public)
        P0; % reference pressure
    end
    
    % The constructor
    methods
        function obj = ConstPressure(P)
            if nargin > 0
                assert(isnumeric(P) && isscalar(P))
                obj.P0 = P;
            end
        end
    end
    
    % Required barotrope methods
    methods
        function test(obj)
            disp(obj)
        end
        
        function P = pressure(obj,~)
            P = obj.P0;
        end
        
        function rho = density(~,~)
            rho = NaN;
        end
    end
    
end

