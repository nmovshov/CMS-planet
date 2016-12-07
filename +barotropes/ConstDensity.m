classdef ConstDensity < barotropes.Barotrope
    %CONSTDENSITY A constant density barotrope.
    
    properties (Access = public)
        rho0; % reference density
    end
    
    % The constructor
    methods
        function obj = ConstDensity(rho)
            if nargin > 0
                assert(isnumeric(rho) && isscalar(rho))
                obj.rho0 = rho;
            end
        end
    end
    
    % Required barotrope methods
    methods
        function PF = test(~)
            PF = true;
        end
        
        function P = pressure(~,rho)
            P = NaN(size(rho));
        end
        
        function rho = density(obj,P)
            rho = obj.rho0*ones(size(P));
        end
    end
    
end

