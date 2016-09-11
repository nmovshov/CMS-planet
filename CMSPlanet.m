classdef CMSPlanet < handle
    %CMSPLANET Interior model of rotating fluid planet.
    %   This class implements a model of a rotating fluid planet using the
    %   Concentric Maclaurin Spheroids technique. A CMSPlanet object is defined
    %   by a given mass, equatorial radius, rotation period, and barotrope.
    
    %% Properties
    properties (Access = public)
        name
        desc
        cms
        baro
        opts
    end
    
    %% The constructor
    methods
        function obj = CMSPlanet(varargin)
            %CMSPLANET Class constructor.
        end
    end
    
end

