%% CONSTANT DENSITY PLANET
% Example and test of the CMSPlanet class. We construct and converge a model of
% a rotating fluid planet with constant density. The cms function handles the
% solution of the shape functions and resulting gravity moments (help cms.m)
% but the preferred workflow when working with physical bodies is to use the
% CMSPlanet class. This interface class will translate between the physical
% parameters defining the planet of interest and the non-dimensional variables
% used in CMS theory.

%% Prepare workspace
clear
clc
close all
earth_mass = 5.9722e24;
G = 6.6743e-11;

%% Construct the exact analytic solution first, because m=f(q) is hard
m = 0.1;  % fairly fast rotation!
lfun = @(x)(3./(2*x.^3)).*((3 + x.^2).*atan(x) - 3*x) - m;
el = fzero(lfun,0.5); % ellipse parameter
b_exact = sqrt(1/(1 + el^2)); % polar radius
s3 = b_exact; % *mean* radius, s^3=b*a^2 (but a=1) we will use this later
xi_exact = @(mu)1./sqrt((1 + (el^2).*(mu.^2)));

%% Take a quick look for sanity check
theta = linspace(0,2*pi);
polax = polarplot(theta, xi_exact(cos(theta)));
polax.DisplayName = 'Maclaurin''s solution';
polax.Parent.ThetaZeroLocation = 'top';
polax.Parent.ThetaDir = 'clockwise';
polax.Parent.ThetaAxisUnits = 'rad';
hold

%% Set up a CMS Planet with constant density and relax to HE
% A CMSPlanet object is defined by a densiy profile rho(a), supplied by the
% user and stored in the column vectors obj.ai and obj.rhoi, indexed from the
% surface in. To complete the definition the user must also specify a mass,
% equatorial radius, and rotation period. With these a gravity field and
% equilibrium shape can be determined, with a call to obj.relax_to_HE().

N = 2; % number of layers (doesn't matter for const. density!)
cmp = CMSPlanet; % an empty object
cmp.ai = linspace(1, 1/N, N); % normalized layer radii
cmp.rhoi = ones(size(cmp.ai)); % normalized layer densities (here constant)
cmp.radius = 7e4*1e3; % ref. equatorial radius of planet [m]
cmp.mass = 317*5.97e24; % ref. mass of real planet [kg]
cmp.renormalize(); % adjust the density and radius to real units
cmp.period = 2*pi/sqrt(m/s3*G*cmp.M/cmp.a0^3); % m=w^2s^3/GM of exact solution above
cmp.opts.verbosity = 2; % usually quieter is better but you have the option...
cmp.relax_to_HE(); % fast for 2 layers; can take a while for high-res models
cmp.renormalize(); % no change to Js!

%% Take a quick look for sanity check
mu = [0, cmp.CMS.mus];
xi_cms = [1, cmp.CMS.zetas(1,:)];
theta = acos(mu);
polarplot(theta, xi_cms, 'DisplayName', 'CMS solution');
legend

%% Compare numerical and analytic solutions (still relies on underlying cms)
% Compare the level surface radii
dxi = xi_cms - xi_exact(mu);
xi_err = max(abs(dxi));
figure
lh = semilogy(mu, abs(dxi));
ah = lh.Parent;
ah.XLabel.Interpreter = 'latex';
ah.YLabel.Interpreter = 'latex';
ah.Title.Interpreter = 'latex';
ah.XLabel.String = '$\mu = \cos(\theta)$';
ah.YLabel.String = '$d\xi$';
ah.Title.String = '$d\xi = \xi(\mu) - 1/\sqrt{1 + l^2\mu^2}$';
fprintf('CMS shape error = %g\n', xi_err)

% Compare the J values
n = 0:2:10;
Js_exact = (-1).^(1 + n/2).*(3./((n + 1).*(n + 3))).*(el^2/(1 + el^2)).^(n/2);
Js = cmp.Js(1:6);
dJ = Js - Js_exact;
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
fprintf('CMS J2 relative error = %g\n',abs(Js(2) - Js_exact(2))/Js_exact(2));

%% The CMSPlanet class has some convenience utilities, for plotting etc.
cmp.plot_rho_of_r();
cmp.report_card()
