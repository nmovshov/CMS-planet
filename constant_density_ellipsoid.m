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

%% Get surface layer shape
a = cms.zeta_j_of_mu(1,0); % should be 1.0!
b = cms.zeta_j_of_mu(1,1);
mu = [0, cms.mus, 1];
xi_cms = [a, cms.zetas(1,:), b];
assert(abs(a - 1) < 1e-12)

%% Construct Maclaurin ellipsoid
% Iterate to converge on Maclaurin ell parameter (?!?!)
q = opts.qrot;
m = 0.001;
for k=1:20
    lfun = @(x)(3./(2*x.^3)).*((3 + x.^2).*atan(x) - 3*x) - m;
    el = fzero(lfun,0.5);
    b_Mac = sqrt(1/(1 + el^2));
    s_Mac = b_Mac; % s=b*a^2, a=1
    m = q*s_Mac^3; % m=q*s^3/a^3, a=1
end
fprintf('m=%g; l=%g\n',m,el)
xi_ell = 1./sqrt((1 + (el.^2).*(mu.^2)));

%% Compare
dxi = xi_cms - xi_ell;
semilogy(mu, abs(dxi))
