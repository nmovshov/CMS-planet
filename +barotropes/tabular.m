classdef tabular < barotropes.Barotrope
    %TABULAR A base class for a generic table barotrope.
    
    %% Properties
    properties (SetAccess = protected)
        mP   % a vector of pressure values
        mrho % a vector of density values
    end
    properties (Access = public)
        interp_method
        extrap_method
        meta
    end
    
    %% The constructor
    methods
        function obj = tabular(P, rho, intm, extm)
            if (nargin == 0) && (nargout == 0)
                print_usage()
                return
            end
            try
                narginchk(2,4)
                assert(isnumeric(P) && isvector(P) && all(double(P) >= 0))
                assert(isnumeric(rho) && isvector(rho) && all(double(rho) >= 0))
                assert(length(P) == length(rho))                
            catch ME
                print_usage()
                rethrow(ME)
            end
            obj.mP = P;
            obj.mrho = rho;
%             if nargin > 2, obj.interp_method = intm; end
%             if nargin > 3, obj.extrap_method = extm; end
%             try
%                 interp1(P, rho, P(1), obj.interp_method, obj.extrap_method);
%             catch ME
%                 print_usage()
%                 rethrow(ME)
%             end
        end
    end
    
    %% Required barotrope methods
    methods
        function test(obj)
            disp(obj)
        end
        
        function P = pressure(obj,rho)
            P = interp1(obj.mrho, obj.mP, rho);
        end
        
        function rho = density(obj,P)
            rho = interp1(obj.mP, obj.mrho, P);
        end
    end
end

%% Usage message
function print_usage()
fprintf('Usage: barotropes.tabular(P, rho)\n')
fprintf('positional arguments:\n')
fprintf('  P    a vector of pressure values (real nonnegative)\n')
fprintf('  rho  a vector of density values(real nonnegative)\n')
end
