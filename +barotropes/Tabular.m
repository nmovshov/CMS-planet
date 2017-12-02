classdef Tabular < barotropes.Barotrope
    %TABULAR A base class for a generic table barotrope.
    
    %% Properties
    properties (SetAccess = protected)
        P_vals   % a vector of pressure values
        rho_vals % a vector of density values
    end
    properties (Access = public)
        interpolation_method
        extrapolation_method
        meta
    end
    
    %% The constructor
    methods
        function obj = Tabular(P, rho, intm, extm)
            if (nargin == 0)
                if nargout > 0, return, end
                print_usage()
                clear obj
                return
            end
            try
                narginchk(2,4)
                if nargin < 4, extm = NaN; end
                if nargin < 3, intm = 'linear'; end
                assert(isnumeric(P) && isvector(P) && all(double(P) >= 0))
                assert(isnumeric(rho) && isvector(rho) && all(double(rho) >= 0))
                assert(length(P) == length(rho),...
                    'len(P) = %i ~= %i = len(rho)', length(P), length(rho))
                assert(interp1(P, rho, P(1), intm, extm) == rho(1))
            catch ME
                print_usage()
                rethrow(ME)
            end
            obj.P_vals = P;
            obj.rho_vals = rho;
            obj.interpolation_method = intm;
            obj.extrapolation_method = extm;
        end
    end
    
    %% Required barotrope methods
    methods
        function PF = test(~)
            PF = true;
        end
        
        function P = pressure(obj,rho)
            P = interp1(obj.rho_vals, obj.P_vals, rho);
        end
        
        function rho = density(obj,P)
            rho = interp1(obj.P_vals, obj.rho_vals, P);
        end
    end
end

%% Usage message
function print_usage()
fprintf('Usage: barotropes.tabular(P, rho, intm, extm)\n')
fprintf('positional arguments:\n')
fprintf('  P    a vector of pressure values [ real nonnegative ]\n')
fprintf('  rho  a vector of density values [ real nonnegative ]\n')
fprintf('  intm (optional) interpolation method [ {''linear''} | ''nearest'' | ''spline'' | ''pchip'' ]\n')
fprintf('  extm (optional) extrapolation strategy [ ''extrap'' | scalar value | {NaN} ]\n')
end
