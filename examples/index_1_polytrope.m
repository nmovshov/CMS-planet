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
%addpath(fullfile(pwd,'..'))
try
    si = setUnits; % if you have physunits
catch
    si = setFUnits; % if you don't have physunits
end
G = si.gravity;

%% Set up a CMS Planet with arbitrary mass and radius
% I'm using numbers for Jupiter just for kicks, no effect on Js of course.
M = 317.8*si.earth_mass;
R = 71492*si.km;

cmp = CMSPlanet(512);
cmp.name = [int2str(cmp.nlayers),'-layer CMS'];
cmp.desc = 'An index-1 polytropic planet';
cmp.M = M;
cmp.a0 = R;
cmp.qrot = 0.089195487; % Hubbard (2013) Table 5

%% Construct a polytrope of index 1 to represent the planet's eos
n = 1;
K = 2*G/pi*R^2; % ...matches radius just for show, K has no effect on the Js
eos = barotropes.Polytrope(K, n);
eos.name = '$P\propto\rho^2$';
cmp.eos = eos;

%% To speed up convergence start with an approximate density structure
a = sqrt(2*pi*G/K);
r = pi/a;
rho_av = 3*M/(4*pi*r^3);
rho_c = (pi^2/3)*rho_av;
x = (cmp.ai(1:end-1) + cmp.ai(2:end))/2;
x(end+1) = cmp.ai(end)/2;
cmp.rhoi = rho_c*sin(a*x)./(a*x);

%% Relax to desired barotrope
cmp.opts.verbosity = 3;
cmp.opts.MaxIterHE = 12;
cmp.opts.dBtol = 1e-12;
cmp.opts.email = '';
cmp.relax_to_barotrope;

%% Compare computed and analytic density structure
q = cmp.qrot;
% Zharkov & Trubistyn (1978) eq. 34.12
ZT3 = [q;...
    (0.173273*q - 0.197027*q^2 + 0.15*q^3)*1e2;...
    (-0.081092*q^2 + 0.15*q^3)*-1e4;...
    (0.056329*q^3)*1e5;...
    nan; nan; nan];
% Hubbard (2013) Table 5
H13_256 = [q; 1.3991574; 5.3182810; 3.0118323; 2.1321157; 1.7406710; 1.5682179];
H13_512 = [q; 1.3989253; 5.3187997; 3.0122356; 2.1324628; 1.7409925; 1.5685327];

% CMSPlanet
CMP = [q; cmp.Js(2)*1e2; -cmp.Js(3)*1e4; cmp.Js(4)*1e5; -cmp.Js(5)*1e6;...
    cmp.Js(6)*1e7; -cmp.Js(7)*1e8;];

% Make it a table
T = table(ZT3, H13_256, H13_512, CMP);
T.Properties.RowNames = {'q','J2x10^2','-J4x10^4','J6x10^5','-J8x10^6',...
    'J10x10^7','-J12x10^8'};

% Display
format long
format compact
disp(T)
format
try
    cmp.plot_equipotential_surfaces;
    cmp.plot_barotrope;
catch
end

%% Save and deliver
% save('index1polytrope', 'cmp', 'T')
% try
% sendmail('address','n=1 polytrope','','index1polytrope.mat')
% catch
% end
% !shutdown /s /t 30
% exit force
