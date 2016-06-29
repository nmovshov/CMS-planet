classdef IdealGasIsotherm < barotropes.Barotrope
    %IDEALGASISOTHERM A toy barotrope - ideal gas at some constant T.
    
    %% Properties
    properties (Access = public)
        T0 % reference temperature (K)
        MMW % mean molecular weight (dimensionless, but numerically = to g/mol)
        molar_mass % kg/mole
    end
    
    properties (Constant)
        R = 8.3145; % J mol^-1 K^-1
    end
    
    %% The constructor
    methods
        function obj = IdealGasIsotherm(mmw, T)
            if nargin > 0
                assert(isnumeric(mmw) && isscalar(mmw))
                assert(isnumeric(T) && isscalar(T))
                obj.MMW = mmw;
                obj.molar_mass = mmw/1000;
                obj.T0 = T;
            end
        end
    end
    
    %% Required barotrope methods
    methods
        function test(obj)
            disp(obj)
        end
        
        function P = pressure(obj,rho)
            n = rho/obj.molar_mass;
            P = n*obj.R*obj.T0;
        end
        
        function rho = density(obj,P)
            n = P/obj.R/obj.T0;
            rho = n*obj.molar_mass;
        end
    end
    
end

