%% CONSTANT DENSITY PLANET
% Example and test of the CMSPlanet class. We will construct and converge a
% model of a rotating fluid planet with constant density. The underlying
% concentric Maclauring spheroids model handles the solution of the shape
% functions and resulting gravity moments (see constant_density_ellipsoid.m) but
% the preffered way to initiallize the CMS object when working with physical
% bodies is to use the CMSPlanet class. This interface class will translate
% between the physical parameters defining the planet of interest and the
% nondimensional variables used in CMS theory.

%% Prepare workspace
clear
clc
close all
addpath(fullfile(pwd,'..'))
try
    si = setUnits; % if you have physunits
catch
    si = setFUnits; % if you don't have physunits
end

%% Set up a CMS Planet with constant density
cmp = CMSPlanet(4); % nlayers shouldn't matter
cmp.a0 = 7e4*si.km;
cmp.M = 317*si.earth_mass;
cmp
