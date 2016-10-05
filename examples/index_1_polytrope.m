%% INTERIOR STRUCTURE OF ROTATING PLANET WITH INDEX 1 POLYTROPE DENSITY
% Example and test of the CMSPlanet class. We construct and converge a
% model of a rotating fluid planet with the pressure-density law
%
% $$P = K\rho^2$$
%
% with a polytropic constant $K$ chosen to match Jupiter's mass and radius.

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

%% Set up a CMS Planet with a given mass and rotation period
% http://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
M = 317.8*si.earth_mass;
R = 71492*si.km;
T = 9.9250*si.hr;
w = 2*pi/T;
G = si.gravity;
q = w^2*R^3/G/M;

cmp = CMSPlanet(4);
cmp.M = M;
cmp.a0 = R;
cmp.qrot = q;

%% Construct a polytrope of index 1 to represent the planet's eos
eos = barotropes.Polytrope(1e5, 1);
cmp.eos = eos;
