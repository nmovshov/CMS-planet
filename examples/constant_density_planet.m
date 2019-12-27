%% CONSTANT DENSITY PLANET
% Example and test of the CMSPlanet class. We will construct and converge a model
% of a rotating fluid planet with constant density. The cms function handles the
% solution of the shape functions and resulting gravity moments (help cms.m) but
% the preferred workflow when working with physical bodies is to use the CMSPlanet
% class. This interface class will translate between the physical parameters
% defining the planet of interest and the non-dimensional variables used in CMS
% theory.

%% Prepare workspace
clear
clc
close all
si = setFUnits;

%% Construct the exact analytic solution first, because m=f(q) is hard
m = 0.1; % small rotation parameter
lfun = @(x)(3./(2*x.^3)).*((3 + x.^2).*atan(x) - 3*x) - m;
el = fzero(lfun,0.5); % ellipse parameter
b_exact = sqrt(1/(1 + el^2)); % polar radius
s3 = b_exact; % *mean* radius, s^3=b*a^2 (but a=1) we will use this later
xi_exact = @(mu)1./sqrt((1 + (el^2).*(mu.^2)));

% Take a quick look for sanity check
fh = figure;
theta = linspace(0,2*pi);
polh = polarplot(theta, xi_exact(cos(theta)));
polh.DisplayName = 'Exact ellipsoid';
polh.Parent.ThetaZeroLocation = 'top';
polh.Parent.ThetaDir = 'clockwise';
polh.Parent.ThetaAxisUnits = 'rad';
polh.DisplayName = 'Exact solution';
hold(polh.Parent, 'on')

%% Set up a CMS Planet with constant density and relax to HE
N = 2; % number of layers doesn't matter for const. density!
cmp = CMSPlanet; % an empty object
cmp.ai = linspace(1, 1/N, N); % normalized layer radii
cmp.rhoi = ones(size(cmp.ai)); % normalized layer densities
cmp.radius = 7e4*si.km; % ref equatorial radius of real planet
cmp.mass = 317*si.earth_mass; % ref mass of real planet
cmp.qrot = m/s3; % CMS method uses q=w^2a^3/GM as rotation parameter
cmp.opts.verbosity = 2; % Usually quieter is better but you have the option...
cmp.relax_to_HE; % This is fast for 2 layers but can take a while for real cases
cmp.renormalize; % Adjust the density and radius to real units, no change to Js!

%% The converged cmp holds the shape information returned by cms.m in field CMS
mu = [0, cmp.CMS.mus];
xi_cms = [1, cmp.CMS.zetas(1,:)];

% Take a quick look for sanity check
polax = polarplot(acos(mu), xi_cms);
polax.DisplayName = 'CMS solution';
polax.Parent.ThetaZeroLocation = 'top';
polax.Parent.ThetaDir = 'clockwise';
polax.Parent.ThetaAxisUnits = 'rad';
legend('show')

%% Compare numerical and analytic solutions (still relies on underlying cms)
figure
% Compare the level surface radii
dxi = xi_cms - xi_exact(mu);
lh = semilogy(mu, abs(dxi));
ah = lh.Parent;
ah.XLabel.Interpreter = 'latex';
ah.YLabel.Interpreter = 'latex';
ah.Title.Interpreter = 'latex';
ah.XLabel.String = '$\mu = \cos(\theta)$';
ah.YLabel.String = '$d\xi$';
ah.Title.String = '$d\xi = \xi(\mu) - 1/\sqrt{1 + l^2\mu^2}$';
xi_err = max(abs(dxi));

% Compare the J values
k = 0:6;
n = 2*k;
J_exact = (-1).^(1 + n/2).*(3./((n + 1).*(n + 3))).*(el^2/(1 + el^2)).^(n/2);
J_cms = cmp.Js(k+1);
dJ = (J_cms - J_exact)./J_exact;
subplot(2,1,1,ah);
subplot(2,1,2);
lh = stem(n, abs(dJ));
ah = lh.Parent;
ah.YScale = 'log';
ah.XTick = n(1:2:end);
ah.XMinorTick = 'on';
ah.XLabel.Interpreter = 'latex';
ah.YLabel.Interpreter = 'latex';
ah.Title.Interpreter = 'latex';
ah.XLabel.String = '$n = 0,2,4,\ldots$';
ah.YLabel.String = '$dJ_n$';
sform = '$dJ_n = J_n - \frac{3(-1)^{1+n/2}}{(n+1)(n+3)}\left(\frac{l^2}{1+l^2}\right)^{n/2}$';
ah.Title.String = sform;
J_err = max(abs(dJ)) %#ok<NOPTS>

%% The CMSPlanet class has some convenience utilities, for plotting etc.
cmp.plot_rho_of_r;
cmp.report_card()
