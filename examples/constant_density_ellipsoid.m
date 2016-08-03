%% MACLAURIN ELLIPSOID WITH CMS
% The Maclaurin spheroid is a closed analytic solution for the shape of a
% rotating *constant density* self-gravitating fluid. This script compares the
% numerical CMS solution with the expected analytic solution.

%% Maclaurin's solution
% The equillibrium shape of a constant density rotating fluid can be shown
% (although not by me) to be an ellipsoid of revolution such that the radius $r$
% follows
% 
% $$r^2(\mu) = \frac{a^2}{1 + l^2\mu^2}$$
% 
% where $a$ is the equatorial radius, $\mu$ is the cosine of the angle from the
% rotation axis (the colatitude), $b$ is the polar radius and
%
% $$l^2 = \frac{a^2}{b^2} - 1.$$
%
% The ellipticity parameter $l$ is related to a dimensionless rotation
% parameter $m=\omega^2s^3/GM$ by the trancendental equation
%
% $$ m = \frac{3}{2l^3}[(3 + l^2)\arctan{l} - 3l]. $$
% 
% The rotation parameter $m$ is given in terms of the ellipsoid's *mean* radius
% $s$. For an oblate ellipsoid the radii are related by $s^3=ba^2$.

%% Prepare workspace
clear
clc
close all
addpath(fullfile(pwd,'..'))

%% Construct an exact normalized (a=1) Maclaurin ellipsoid
m = 0.1; % small rotation parameter
lfun = @(x)(3./(2*x.^3)).*((3 + x.^2).*atan(x) - 3*x) - m;
el = fzero(lfun,0.5); % ellipse parameter
b = sqrt(1/(1 + el^2)); % polar radius
s3 = b; % *mean* radius, s^3=b*a^2 (but a=1) we will use this later
xi_exact = @(mu)1./sqrt((1 + (el^2).*(mu.^2)));

% Let's take a quick look
theta = linspace(0,2*pi);
polax = polarplot(theta, xi_exact(cos(theta)));
polax.Parent.ThetaZeroLocation = 'top';
polax.Parent.ThetaDir = 'clockwise';
polax.Parent.ThetaAxisUnits = 'rad';


%% Set up a CMS object to mimic constant density case
opts = cmsset;
opts.nlayers = 1;
opts.rcore = 0.01;
opts.qrot = 0.1;
opts.kmax = 32;
opts.dJtol = 1e-14;
opts.verbosity = 2;
cms = ConcentricMaclaurinSpheroids(opts);

%% Converge cms (this may take a while)
cms.relax

%% Get surface layer shape
a = cms.zeta_j_of_mu(1,0); % should be 1.0!
b = cms.zeta_j_of_mu(1,1);
mu = [0, cms.mus, 1];
xi_cms = [a, cms.zetas(1,:), b];
assert(abs(a - 1) < 1e-12)

%% Compare
dxi = xi_cms - xi_exact;
semilogy(mu, abs(dxi))
