classdef ConstPressure < barotrope.Barotrope
    %CONSTPRESSURE A test implementation of constant pressure barotrope.
    
    properties
        P0 = 0;
    end
    
    % The constructor
    methods
        function obj = ConstPressure(p)
            if nargin > 0
                assert(isscalar(p))
                assert(isnumeric(p))
                obj.P0 = p;
            end
        end
    end
    
    % Required barotrope methods
    methods
        function test(obj)
            disp(obj)
        end
        
        function y = P(obj,~)
            y = obj.P0;
        end
        
        function y = rho(~,~)
            y = NaN;
        end
    end
    
end

