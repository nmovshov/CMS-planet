%% CONSTANT DENSITY PLANET
% Example and test of the CMSPlanet class. We will construct and converge a
% model of a rotating fluid planet with constant density. The underlying
% concentric Maclauring spheroids model handles the solution of the shape
% functions and resulting gravity moments (see constant_density_ellipsoid.m) but
% the preffered way to initiallize the CMS object when working with physical
% bodies is to use the CMSPlanet class. This interface class will translate
% between the physical parameters defining the planet of interest and the
% nondimensional variables used in CMS theory.

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

%% Construct the exact analytic solution first, because m=f(q) is hard
m = 0.1; % small rotation parameter
lfun = @(x)(3./(2*x.^3)).*((3 + x.^2).*atan(x) - 3*x) - m;
el = fzero(lfun,0.5); % ellipse parameter
b_exact = sqrt(1/(1 + el^2)); % polar radius
s3 = b_exact; % *mean* radius, s^3=b*a^2 (but a=1) we will use this later
xi_exact = @(mu)1./sqrt((1 + (el^2).*(mu.^2)));

% Take a quick look for sanity check
try % requires R2016a or later
    fh = figure;
    theta = linspace(0,2*pi);
    polh = polarplot(theta, xi_exact(cos(theta)));
    polh.Parent.ThetaZeroLocation = 'top';
    polh.Parent.ThetaDir = 'clockwise';
    polh.Parent.ThetaAxisUnits = 'rad';
    polh.DisplayName = 'Exact solution';
catch
    delete(fh);
end

%% Set up a CMS Planet with constant density and relax to HE
cmp = CMSPlanet(4); % nlayers shouldn't matter
cmp.a0 = 7e4*si.km;
cmp.M = 317*si.earth_mass;
cmp.rhoi = ones(cmp.nlayers,1)*cmp.rho0; % for illustration
cmp.qrot = m/s3; % CMS method uses q=w^2a^3/GM as rotation parameter
cmp.cms.relax();

%% Compare numerical and analytic solutions (still relies on underlying cms)

% Compare the level surface radii
b_cms = cmp.cms.bs(1);
mu = [0, cmp.cms.mus, 1];
xi_cms = [1, cmp.cms.zetas(1,:), b_cms];
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
b_err = b_cms - b_exact %#ok<NOPTS>

% Compare the J values
n = 0:2:cmp.opts.kmax;
J_exact = (-1).^(1 + n/2).*(3./((n + 1).*(n + 3))).*(el^2/(1 + el^2)).^(n/2);
J_cms = cmp.Jn;
dJ = J_cms - J_exact;
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
