%% GRAVITY COEFFICIENTS OF ROTATING INDEX-1 POLYTROPE
% Example and test of the CMSPlanet class. We construct and converge a
% model of a rotating fluid planet with the pressure-density law
%
% $$P = K\rho^2$$
%
% with a polytropic constant $K$. This script demonstrates how to set up a
% CMSPlanet object with a specified barotrope and converge to a density
% structure in hydrostatic equilibrium with the given barotrope. The default
% starting density profile is that of a homogeneous sphere.
%
% The comparison model is the 3rd-order Zharkov and Trubitsyn theory. I take the
% numbers from Table 5 of Hubbard (2013).

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

cmp = CMSPlanet(16);
cmp.name = 'poly1';
cmp.desc = 'An index-1 polytropic planet';
cmp.M = M;
cmp.a0 = R;
cmp.qrot = 0.089195487; % used in Hubbard (2013)

%% Construct a polytrope of index 1 to represent the planet's eos
% The initial value of K is chosen to match the given radius if the planet were
% not rotating. We cannot expect an exact match when we specify mass, radius,
% and eos at the same time as this over determines the planet's structure. The
% CMS method will try to match the density structure to the barotrope and
% equatorial radius, leaving the total mass different from the specified value.
% Then the polytropic constant is adjusted to match the total mass. I _think_
% this is what is done in Hubbard (2013).
n = 1; % polytrope index
K = 2*G/pi*R^2; % polytrope constant
eos = barotropes.Polytrope(K, n);
cmp.eos = eos;

%% Relax to desired barotrope
cmp.opts.verbosity = 1;
cmp.opts.dBtol = 0.0001;
cmp.opts.MaxIterBar = 10;
cmp.relax_to_barotrope;

%% Compare computed and analytic density structure
q = cmp.qrot;
ZT3 = [q; 1.3994099; 5.3871087; 3.9972442; nan; nan; nan];
H13_256 = [q; 1.3991574; 5.3182810; 3.0118323; 2.1321157; 1.7406710; 1.5682179];
H13_512 = [q; 1.3989253; 5.3187997; 3.0122356; 2.1324628; 1.7409925; 1.5685327];
CMP = [q; cmp.Js(2)*1e2; -cmp.Js(3)*1e4; cmp.Js(4)*1e5; -cmp.Js(5)*1e6; cmp.Js(6)*1e7; -cmp.Js(7)*1e8;];
T = table(ZT3, H13_256, H13_512, CMP);
T.Properties.RowNames = {'q','J2x10^2','-J4x10^4','J6x10^5','-J8x10^6','J10x10^7','-J12x10^8'};
