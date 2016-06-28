classdef (Abstract) Barotrope
    %BAROTROPE Interface base class for all P(rho) relations.
    %   The BAROTROPE abstract base class will implement only initialization of
    %   named constants and define syntax for interface to specialized
    %   subclasses.
    
    properties
        si;
    end
    
    % Constructor
    methods
        function obj = Barotrope()
            % This constructor is called implicitly with no arguments when a
            % subclass instance is created.
        end
    end
    
    % Required methods for any derived subclass
    methods (Abstract)
        test(obj)
    end
end
