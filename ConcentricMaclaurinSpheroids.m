classdef ConcentricMaclaurinSpheroids
    %CONCENTRICMACLAURINSPHEROIDS Implementation of CMS shape model.
    %   This class implements the iterative relaxation of concentric Maclaurin
    %   spheroids with given (dimensionless) radii and densities to a
    %   self-consistent hydrostatic equilibrium shape, as explained in Hubbard
    %   (2013).
    
    %% Properties
    properties (Access = public)
        opts
    end
    properties (GetAccess = public, SetAccess = private)
        lambdas
        deltas
        zetas
        Js
    end
    
    %% The constructor
    methods
        function obj = ConcentricMaclaurinSpheroids(opts)
            if nargin == 0
                opts = ConcentricMaclaurinSpheroids.cmsset();
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
        function options = cmsset(varargin)
        %CMSSET Create cms options structure.
        %   OPTIONS = CMSSET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an
        %   options structure OPTIONS in which the named properties have the
        %   specified values. Any unspecified properties have default values.
        %   Case is ignored for property names.
        %
        %   CMSSET with no input arguments displays all property names and their
        %   possible values.
        %
        %CMSSET PROPERTIES
        %
        %nlayers - Number of constant density layers [positive integer {512}]
        %
        %nangles - Number of colatitude points used to define level surfaces [positive integer {48}]
        %
        %nmoments - Degree to carry out mulitpole expansion of gravity moments [positive even {12}]
        %
        %   Note: defaults chosen to match Hubbard (2013) example.
        
            % Print out possible values of properties.
            if (nargin == 0) && (nargout == 0)
                fprintf('  nlayers: [positive scalar integer {512}]\n');
                fprintf('  nangles: [positive scalar integer {48}]\n');
                fprintf('  nmoments: [positive scalar even {12}]\n');
                return
            end
            
            % Parse name-value pairs.
            if rem(nargin, 2) ~= 0
                error('Name/Value pair mismatch')
            end
            %TODO:CONTINUE
            
        end
        
    end % End of static methods
    
end % End of classdef

