%% A PLANET WITH POLYTROPIC ENVELOPE AND CONSTANT DENSITY CORE
% This example can be used as a template to construct an interior model of a
% planet based on a polytrope with a constant density core.

%% Prepare workspace
clear
clc
close all
%addpath(fullfile(pwd,'..'))
try
    si = setUnits; % <FEX>/13018-physunits-module-from-fortran
catch
    si = setFUnits; % if you don't have physunits
end
G = si.gravity;

%% Specify required parametrs (if using as template modify this section only)
% Mass and radius
M = 317.8*si.earth_mass; % total planet mass
R = 71492*si.km;         % equatorial radius

% Rotation paramater
P_rot = 9*si.hour + 55*si.minute + 29.7*si.second;
w_rot = 2*pi/P_rot;
q = w_rot^2*R^3/(G*M); % (or you can specify q directly)

% Polytrope index n and constant K
poly_n = 1;
poly_K = 2e5*si.Pa/(si.kg/si.m^3)^(1 + 1/poly_n);

% Core properties
r_core = 0.15*R;          % core equatorial radius
m_core = 5*si.earth_mass; % core mass

% Model options
N = 12;          % number of spheroids to use, including core
name = '';       % optional name
desc = '';       % optional description
opts = cmsset(); % use for all other model options, see help cmsset for a list

%% Create the CMSPlanet object
cmp = CMSPlanet(N, opts);
cmp.name = name;
cmp.desc = desc;
cmp.M = M;
cmp.a0 = R;
cmp.qrot = q;

%% Set layer radii using N-1 layers for envelope and one for core
cmp.ai = linspace(1, r_core/cmp.a0, cmp.nlayers)*cmp.a0;

%% Assign an array of (handles to) barotropes to the CMSPlanet
eos_env = barotropes.Polytrope(poly_K, poly_n);
eos_env.name = ['$P\propto\rho^',num2str(1 + 1/poly_n),'$'];
eos_core = barotropes.ConstDensity(m_core/(4*pi/3*r_core^3));
eos_core.name = '$\rho=rho_c$';
cmp.eos = [repmat(eos_env,N-1,1); eos_core];

%% To speed up convergence try to start with an approximate density structure
% We use the density profile of a non-rotating n=1 polytrope but we have to
% guard against negative densities and mismatched dimensions. The benefit is a
% skipping a few iterations of relaxation loop. Comment out this section to
% start from a uniform density structure.
a = double(sqrt(2*pi*G/poly_K));
r = pi/a;
rho_av = double(3*M/(4*pi*r^3));
rho_c = (pi^2/3)*rho_av;
x = double(cmp.ai(1:end-1) + cmp.ai(2:end))/2;
x(end+1) = double(cmp.ai(end))/2;
rho_guess = rho_c*(sin(a*x)./(a*x));
rho_guess(rho_guess < 0) = 0;
cmp.rhoi = rho_guess*cmp.rho0;
cmp.rhoi = cmp.rhoi*cmp.M/cmp.M_calc;

%% Relax to desired barotrope
cmp.opts.email = '';
cmp.relax_to_barotrope;

%% Display
try
    cmp.plot_equipotential_surfaces;
    cmp.plot_barotrope('showinput');
catch
end

%% (Optional) Save and deliver
% save('filename', 'cmp')
% try
% sendmail('address','subject','msg','filename.mat')
% catch
% end
% !shutdown /s /t 30
% exit force
