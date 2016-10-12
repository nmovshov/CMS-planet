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

%% Set up a CMS Planet with arbitrary mass and radius
M = 317.8*si.earth_mass;
R = 71492*si.km;

cmp = CMSPlanet(16);
cmp.name = 'nrpoly';
cmp.desc = 'A nonrotating polytropic planet';
cmp.M = M;
cmp.a0 = R;
cmp.qrot = 0;
cmp.opts.verbosity = 0;
cmp.opts.MaxIterBar = 2;

%% Construct a polytrope of index 1 to represent the planet's eos
K = 2e5*si.kg^-1*si.m^5*si.s^-2;
n = 1;
eos = barotropes.Polytrope(K, n);
cmp.eos = eos;
