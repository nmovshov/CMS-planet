classdef (Abstract) Barotrope
    %BAROTROPE Interface base class for all barotropes.
    %   The BAROTROPE abstract base class will implement only initialization of
    %   named constants and define syntax for interface to specialized
    %   subclasses.
    
    properties (Constant)
        
    end
    
    % Constructor
    methods
        function obj = Barotrope()
        % ABC constructor called implicitly with no arguments by subclass.
            
        end
    end
    
    % Required methods for any derived subclass
    methods (Abstract)
        test(obj)
        P(obj,rho)
        rho(obj,P)
    end
end
