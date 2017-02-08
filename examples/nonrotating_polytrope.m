%% INTERIOR STRUCTURE OF A NON-ROTATING INDEX 1 POLYTROPE
% Example and test of the CMSPlanet class. We construct and converge a
% model of a non-rotating fluid planet with the pressure-density law
%
% $$P = K\rho^2$$
%
% with a polytropic constant $K$. This script demonstrates how to set up a
% CMSPlanet object with a specified barotrope and converge to a density
% structure in hydrostatic equilibrium with the given barotrope. The default
% starting density profile is that of a homogeneous sphere. This scripts also
% explores how increasing the number of layers reduces the discretization error.
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
%addpath(fullfile(pwd,'..'))
try
    si = setUnits; % if you have physunits
catch
    si = setFUnits; % if you don't have physunits
end
G = si.gravity;

%% Set up some CMS Planets with arbitrary mass and radius
M = 317.8*si.earth_mass;
R = 71492*si.km;

%#ok<*SAGROW>
%nbl = [16, 64, 256];
%nbl = [8, 16];
nbl = [16, 32, 128];
for k=1:length(nbl)
    cmp(k) = CMSPlanet(nbl(k));
    cmp(k).name = ['CMS',int2str(nbl(k))];
    cmp(k).desc = 'A non-rotating polytropic planet';
    cmp(k).M = M;
    cmp(k).a0 = R;
    cmp(k).qrot = 0;
end

%% Construct a polytrope of index 1 to represent the planet's eos
% A non-rotating, index-1 polytrope is completely defined by K. The radius is
% independent of mass. The density structure normalized to the central density
% is also independent of mass, and the absolute value of density is determined
% with the known average density.
n = 1; % polytrope index
K = 2*G/pi*R^2; % polytrope constant
eos = barotropes.Polytrope(K, n);
eos.name = sprintf('$P\\propto\\rho^2$');
[cmp.eos] = deal(eos);

%% Relax to desired barotrope
for k=1:length(nbl)
    cmp(k).opts.verbosity = 1;
    cmp(k).opts.dBtol = 1e-10; % default 1e-10
    cmp(k).opts.MaxIterBar = 40; % default 40
    cmp(k).relax_to_barotrope;
end

%% Compare computed and analytic density structure
% prepare
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesBox', 'on')

% calculate
a = sqrt(2*pi*G/K);
R = pi/a;
rho_av = 3*M/(4*pi*R^3);
rho_c = (pi^2/3)*rho_av;
r = linspace(0,1)*R;
rho_exact = rho_c*sin(a*r)./(a*r);
rho_exact(1) = rho_c;

% plot
ah = axes; hold(ah, 'on');
l1 = plot(r/R, rho_exact/rho_c, 'k--');
l1.DisplayName = '$\sin(ar)/(ar)$';
for k=1:length(nbl)
    try
        l2(k) = stairs(cmp(k).ai/R, cmp(k).rhoi/cmp(k).betanorm/rho_c, '-');
    catch
        l2(k) = stairs(cmp(k).ai/R, cmp(k).rhoi/rho_c, '-'); % back support
    end
    l2(k).DisplayName = ['$\rho_i/\rho_c$ ',cmp(k).name];
end

% annotate
xlabel('$r/R$')
ylabel('$\rho/\rho_c$')

lh = legend(ah, 'show');
lh.FontSize = 12;

% errors
n = length(nbl);
P_err = (cmp(n).P_c - cmp(n).eos.pressure(rho_c))/cmp(n).eos.pressure(rho_c);
s_tit = sprintf(['$N_\\mathrm{layers}=%g$; ',...
    '$\\beta=%g$; $\\Delta P_c=%g\\%%$'],...
    cmp(n).nlayers, double(cmp(n).betanorm), double(P_err)*100);
title(s_tit)

% polytrope definitions
s_poly{1} = '$a=\sqrt{2\pi{}G/K}$';
s_poly{2} = '$R=\pi/a$';
s_poly{3} = '$\rho_c = \frac{\pi{M}}{4R^3}$';
th = annotation('textbox');
th.Interpreter = 'latex';
th.String = s_poly;
th.FontSize = 12;
th.Position = [lh.Position(1), lh.Position(2)-th.Position(4),...
    lh.Position(3), th.Position(4)];

% Also plot converged and input barotrope
cmp(end).plot_barotrope('showinput',true,'showscaledinput',true);
title(s_tit);
