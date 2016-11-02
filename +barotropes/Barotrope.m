classdef (Abstract) Barotrope < handle
    %BAROTROPE Interface base class for all barotropes.
    %   The BAROTROPE abstract base class will implement only initialization of
    %   named constants and define syntax for interface to specialized
    %   subclasses.
    
    %% Properties
    properties (Constant)
        
    end
    
    %% Constructor
    methods
        function obj = Barotrope()
        % ABC constructor called implicitly with no arguments by subclass.
            
        end
    end
    
    %% Required methods for any derived subclass
    % NOTE: MATLAB does NOT enforce the signature of abstract methods on the
    % corresponding concrete method in subclasses. However it is best practice
    % to pretend that it does.
    methods (Abstract)
        test(obj)
        pressure(obj, rho)
        density(obj, P)
    end
end
