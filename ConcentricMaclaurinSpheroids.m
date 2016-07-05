classdef ConcentricMaclaurinSpheroids
    %CONCENTRICMACLAURINSPHEROIDS Implementation of CMS shape model.
    %   This class implements the iterative relaxation of concentric Maclaurin
    %   spheroids with given (dimensionless) radii and densities to a
    %   self-consistent hydrostatic equilibrium shape, as explained in Hubbard
    %   (2013).
    
    %% Properties
    properties (GetAccess = public, SetAccess = private)
        opts
        lambdas
        deltas
        zetas
        Js
    end
    
    %% The constructor
    methods
        function obj = ConcentricMaclaurinSpheroids(opts)
            if nargin == 0
                opts = cmsset();
            end
            obj.lambdas = linspace(1, 1/opts.nlayers, opts.nlayers)';
            obj.zetas = repmat(obj.lambdas, 1, opts.nangles);
            obj.opts = opts;
        end
    end
    
    %% Ordinary methods
    methods
        %function y = foo(obj, x) ...
    end
    
    %% Static methods
    methods (Static)
        
    end % End of static methods
    
end % End of classdef

