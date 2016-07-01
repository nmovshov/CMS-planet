classdef ConcentricMaclaurinSpheroids
    %CONCENTRICMACLAURINSPHEROIDS Implementation of CMS shape model.
    %   This class implements the iterative relaxation of concentric Maclaurin
    %   spheroids with given (dimensionless) radii and densities to a
    %   self-consistent hydrostatic equilibrium shape, as explained in Hubbard
    %   (2013).
    
    %% Properties
    properties (GetAccess = public, SetAccess = private)
        nlayers
        nangles
        lambda
        delta
        zeta
    end
    
    %% The constructor
    methods
        function obj = ConcentricMaclaurinSpheroids(nlayers, nangles)
            if nargin > 0
                obj.nlayers = nlayers;
                obj.nangles = nangles;
                obj.lambda = linspace(1, 0, nlayers)';
                obj.zeta = repmat(obj.lambda, 1, nangles);
            end
        end
    end
    
end

