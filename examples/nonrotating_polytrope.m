%% INTERIOR STRUCTURE OF A NON-ROTATING INDEX 1 POLYTROPE
% Example and test of the CMSPlanet class. We construct and converge a
% model of a non-rotating fluid planet with the pressure-density law
%
% $$P = K\rho^2$$
%
% with a polytropic constant $K$. This script demonstrates how to set up a
% CMSPlanet object with a specified barotrope and converge to a density
% structure in hydrostatic equilibrium with the given barotrope. The default
% starting density profile is that of a homogeneous sphere.
%
% Reminder: the density as a function of radius of a non-rotating, index-1
% polytropic planet in hydrostatic equilibrium is
% 
% $$\rho(r) = \rho_c \sin(ar)/(ar)$$
%
% where
% 
% $$a = \sqrt{2\pi{G}/K}$$
% 
% and
% 
% $$\rho_c = 3M/(4\pi{R^3})$$
%
% where $M$ is the planet's mass and $R = \pi/a$, which is independent of mass.
% 

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
cmp.name = 'nrpoly';
cmp.desc = 'A nonrotating polytropic planet';
cmp.M = M;
cmp.a0 = R;
cmp.qrot = 0;

%% Construct a polytrope of index 1 to represent the planet's eos
% A non-rotating, index-1 polytrope is completely defined by K. The radius is
% independent of mass. The density structure normalized to the central density
% is also independent of mass, and the absolute value of density is determined
% with the known average density.
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
% prepare
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesBox', 'on')

% calculate
a = sqrt(2*pi*G/K);
R = pi/a;
rho_av = 3*cmp.M_calc/(4*pi*R^3);
rho_c = (pi^2/3)*rho_av;
r = linspace(0,1)*cmp.a0;
rho_exact = rho_c*sin(a*r)./(a*r);
rho_exact(1) = rho_c;

% plot
ah = axes; hold(ah);
l1 = plot(r/cmp.a0, rho_exact/rho_c);
l2 = stairs(cmp.ai/cmp.a0, cmp.rhoi/rho_c, '-');

% annotate
l1.DisplayName = '$\sin(ar)/(ar)$';
l2.DisplayName = '$\rho_i/\rho_c, i=1,\ldots{},N$';

xlabel('$r/a_0$')
ylabel('$\rho/\rho_c$')

legend(ah, 'show')

% errors
P_err = (cmp.P_c - cmp.eos.pressure(rho_c))/cmp.eos.pressure(rho_c);
M_err = (cmp.M_calc - cmp.M)/cmp.M;
s_tit = sprintf(['$N_\\mathrm{layers}=%g$; ',...
    '$\\Delta M_\\mathrm{tot}=%g\\%%$; $\\Delta P_c=%g\\%%$'],...
    cmp.nlayers, double(M_err)*100, double(P_err)*100);
title(s_tit)
