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
        function test(obj)
            disp(obj)
        end
        
        function P = pressure(~,~)
            P = NaN;
        end
        
        function rho = density(obj,~)
            rho = obj.rho0;
        end
    end
    
end

