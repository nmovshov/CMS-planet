classdef ConstPressure < barotropes.Barotrope
    %CONSTPRESSURE A constant pressure barotrope.
    
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
        
        function P = pressure(obj,rho)
            P = obj.P0*ones(size(rho));
        end
        
        function rho = density(~,P)
            rho = NaN(size(P));
        end
    end
    
end

