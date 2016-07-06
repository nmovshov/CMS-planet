classdef ConcentricMaclaurinSpheroids
    %CONCENTRICMACLAURINSPHEROIDS Implementation of CMS shape model.
    %   This class implements the iterative relaxation of concentric Maclaurin
    %   spheroids from starting (dimensionless) radii and densities to a
    %   self-consistent hydrostatic equilibrium shape, as explained in Hubbard
    %   (2013).
    
    %% Properties
    properties (GetAccess = public, SetAccess = private)
        opts % holds CMS project-wide options
        lambdas % normalized layer equatorial radii
        deltas % normalized density steps
        zetas % normalized and scaled level-surface radii
        Js % rescaled dimensionless gravity moments
    end
    
    %% The constructor
    methods
        function obj = ConcentricMaclaurinSpheroids(opts)
            if nargin == 0
                opts = cmsset();
            end
            obj.lambdas = linspace(1, 1/opts.nlayers, opts.nlayers)';
            obj.zetas = ones(opts.nlayers, opts.nangles);
            obj.Js.tilde = zeros(opts.nlayers, opts.nmoments);
            obj.Js.tilde_prime = zeros(opts.nlayers, opts.nmoments);
            obj.Js.pprime = zeros(opts.nlayers, 1);
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

