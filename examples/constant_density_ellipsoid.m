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
b_exact = sqrt(1/(1 + el^2)); % polar radius
s3 = b_exact; % *mean* radius, s^3=b*a^2 (but a=1) we will use this later
xi_exact = @(mu)1./sqrt((1 + (el^2).*(mu.^2)));

% Take a quick look for sanity check
try % requires R2016a or later
    theta = linspace(0,2*pi);
    polax = polarplot(theta, xi_exact(cos(theta)));
    polax.Parent.ThetaZeroLocation = 'top';
    polax.Parent.ThetaDir = 'clockwise';
    polax.Parent.ThetaAxisUnits = 'rad';
catch
end

%% Now set up a CMS object to mimic constant density case
q = m/s3; % CMS method uses q=w^2a^3/GM as rotation parameter
opts = cmsset('qrot', q,...
              'nlayers', 1,... % can be any number though!
              'kmax', 32,...
              'verbosity', 2);
cms = ConcentricMaclaurinSpheroids(opts);
cms.deltas = [1, 0*(2:opts.nlayers)];

%% Converge cms (this may take a minute if using many layers)
cms.relax

%% Get surface layer shape
a_cms = cms.zeta_j_of_mu(1,0); % should be 1.0!
assert(abs(a_cms - 1) < 1e-12)
b_cms = cms.zeta_j_of_mu(1,1);
mu = [0, cms.mus, 1];
xi_cms = [a_cms, cms.zetas(1,:), b_cms];

% Take a quick look for sanity check
try % requires R2016a or later
    polax = polarplot(acos(mu), xi_cms);
    polax.Parent.ThetaZeroLocation = 'top';
    polax.Parent.ThetaDir = 'clockwise';
    polax.Parent.ThetaAxisUnits = 'rad';
catch
end

%% Compare numerical and analytic solutions
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
b_err = b_cms - b_exact;
