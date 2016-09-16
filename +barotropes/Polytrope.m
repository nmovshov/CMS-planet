classdef Polytrope < barotropes.Barotrope
    %POLYTROPE A barotrope of the form P = K*rho^(1 + 1/n).
    
    %% Properties
    properties (Access = public)
        K % polytropic constant, dimensions depend on index
        n % polytropic index, dimensionless
        alpha % 1 + 1/n
    end
    
    %% The constructor
    methods
        function obj = Polytrope(K, n)
            if nargin > 0
                assert(isnumeric(K) && isscalar(K) && double(K) > 0)
                assert(isnumeric(n) && isscalar(n) && double(n) > 0)
                obj.K = K;
                obj.n = n;
                obj.alpha = 1 + 1/n;
            end
        end
    end
    
    %% Required barotrope methods
    methods
        function test(obj)
            disp(obj)
        end
        
        function P = pressure(obj,rho)
            P = obj.K*rho.^(obj.alpha);
        end
        
        function rho = density(obj,P)
            rho = (P/obj.K).^(1/obj.alpha);
        end
    end
    
end

