%% INTERIOR STRUCTURE OF A NONROTATING INDEX 1 POLYTROPE
% Example and test of the CMSPlanet class. We construct and converge a
% model of a nonrotating fluid planet with the pressure-density law
%
% $$P = K\rho^2$$
%
% with a polytropic constant $K$.

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
G = si.gravity;

%% Set up a CMS Planet with arbitrary mass and radius
M = 317.8*si.earth_mass;
R = 71492*si.km;

cmp = CMSPlanet(32);
cmp.name = 'nrpoly';
cmp.desc = 'A nonrotating polytropic planet';
cmp.M = M;
cmp.a0 = R;
cmp.qrot = 0;

%% Construct a polytrope of index 1 to represent the planet's eos
% A non-rotating, index-1 polytrope is completely defined by K. The radius is
% independent of mass. The density structure normalized to the central density
% is also independent of mass, and the absolute valaue of density is determined
% with the known average density.
n = 1; % polytrope index
K = 2*G/pi*R^2; % polytrope constant
eos = barotropes.Polytrope(K, n);
cmp.eos = eos;

%% Relax to desired barotrope
cmp.opts.verbosity = 2;
cmp.opts.dMtol = 0.1;
cmp.opts.MaxIterBar = 10;
cmp.relax_to_barotrope;
