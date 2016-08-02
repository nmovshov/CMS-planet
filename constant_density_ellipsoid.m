%% Maclaurin ellipsoid
% Use CMS to reproduce the constant-density Maclaurin ellipsoid

%% Prepare workspace
clear
clc
close all

%% Set up a CMS object to mimic constant density case
opts = cmsset;
opts.nlayers = 2;
opts.rcore = 0.0001;
opts.qrot = 0.1;
opts.MaxIter = 40;
opts.kmax = 30;
opts.verbosity = 2;
cms = ConcentricMaclaurinSpheroids(opts);

%% Converge cms (this may take a while)
cms.relax

%% Get top level surface
a = cms.zeta_j_of_mu(1,0); % should be 1.0!
b = cms.zeta_j_of_mu(1,1);
mu = [0, cms.mus, 1];
xi_cms = [a, cms.zetas(1,:), b];

%% Construct ellipsoid
el = sqrt((a/b)^2 - 1);
xi_ell = 1./sqrt((1 + (el.^2).*(mu.^2)));

%% Compare
dxi = xi_cms - xi_ell;
semilogy(mu, abs(dxi))
