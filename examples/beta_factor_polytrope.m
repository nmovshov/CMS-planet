%% MASS RESCALING OF A PLANET WITH A POLYTROPE
% Example and test of the CMSPlanet class. We construct and converge a
% model of a rotating fluid planet with the pressure-density law
%
% $$P = K\rho^2$$
%
% with a polytropic constant $K$. This script demonstrates how to rescale
% the polytropic constant $K$ accounting for the mass renormalization
% (beta factor) that we obtain at the end of the process.

%% Prepare workspace
clear
clc
close all
si = setFUnits;
G = si.gravity;

%% Specify required parameters
% Model options
N = 32;          % number of spheroids to use, including core
name = '';       % optional name
desc = '';       % optional description
opts = cmsset(); % use for all other model options, see help cmsset for a list

% Mass and radius
M = 317.8 *si.earth_mass; % total planet mass
R = 71492 *si.km;         % equatorial radius

% Rotation paramater
P_rot = 9 *si.hour + 55 *si.minute + 29.7 *si.second;
w_rot = 2*pi/P_rot;
q = w_rot^2*R^3/(G*M); % (or you can specify q directly)

% Polytrope index n and constant K
poly_n = 1.0;
poly_K = 2*G/pi*R^((3-poly_n)/poly_n)*M^((poly_n-1)/poly_n);

%% Start
accepted = 0;
t_str = tic;
for ic = 1:2;
    if ic == 1;
        fprintf('Starting initial computation...\n\n')
    else
        % Rescale properly the polytrope constant
        fprintf('Rescaling the input polytropic constant by the beta factor...\n\n')
        poly_K = poly_K * beta^(-(1/poly_n));
    end
    
    fprintf('Polytrope K = %.4f\n', poly_K)
    
    % Create the initial CMSPlanet object
    cmp = models.single_polytrope(N, [poly_K, poly_n] );
    cmp.name = name;
    cmp.desc = desc;
    cmp.M = M;
    cmp.a0 = R;
    cmp.qrot = q;
    
    % To speed up convergence try to start with an approximate density structure
    a = pi/R;
    rho_av = double(3*M/(4*pi*R^3));
    rho_c = (pi^2/3)*rho_av;
    x = double(cmp.ai(1:end-1) + cmp.ai(2:end))/2;
    x(end+1) = double(cmp.ai(end))/2; %#ok<SAGROW>
    rho_guess = rho_c*(sin(a*x)./(a*x));
    rho_guess(rho_guess < 0) = min(abs(rho_guess));
    cmp.rhoi = rho_guess*cmp.rho0;
    cmp.rhoi = cmp.rhoi*cmp.M/cmp.M_calc;
    
    % Relax to desired barotrope
    cmp.opts.verbosity = 0;
    cmp.opts.MaxIterBar = 40;
    cmp.opts.dBtol = 1e-8;
    cmp.opts.email = '';
    cmp.opts.equipotential_squeeze = 'mean';
    cmp.relax_to_barotrope;
    
    Js_old = cmp.cms.Js;
    rhoi_old = cmp.rhoi;
    beta = cmp.betanorm;
    
    fprintf('J2 x 10^6 = %.2f \n', cmp.J2*1e6)
    fprintf('Beta factor = %.4f\n\n', beta)
end
