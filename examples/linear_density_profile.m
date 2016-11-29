%% GRAVITY COEFFICIENTS OF LINEAR DENSITY PLANET
% Example and test of the CMSPlanet class. We construct and converge a
% model of a rotating fluid planet with a density profile that is linear in the
% mean radius of a level surface.
%
% The comparison is with a 5th-order Zharkov and Trubitsyn model. I take the
% numbers from Table 1 of Hubbard (2013). A complication is that the Z+T theory
% works in powers of the small parameter m rather than q. A CMS-based model
% would have to iteratively converge on a the desired values of m and s_i. But I
% will use the value of q already found by Hubbard.

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

%% Set up a CMS object and give it a density profile linear in lambda to start
% A constant lambda step (default) calls for constant delta step
cms = ConcentricMaclaurinSpheroids(128);
cms.qrot = 0.088822426; % Hubbard (2013) Table 1
cms.deltas = ones(cms.nlayers, 1);
cms.deltas(1) = 0;

%% Relax to hydrostatic equilibrium
cms.opts.verbosity = 2;
cms.opts.MaxIterHE = 40;
cms.opts.email = '';
cms.relax;

%% After initial relaxation, we iteratively fix deltas ss and re-relax
iter = 0;
for iter=1:20
    fprintf('\n  Fixing density profile to mean radii - iteration %i\n', iter)
    new_deltas = [0; -diff(cms.ss/cms.ss(1))];
    new_deltas = new_deltas/max(new_deltas);
%    new_deltas = [0; cms.deltas(2:end).*diff(cms.ss)./diff(cms.lambdas)]; %H13?
    delta_deltas = sum(abs(new_deltas - cms.deltas));
    fprintf('  deltas array modified by < %g.\n\n',delta_deltas)
    cms.deltas = new_deltas;
    if delta_deltas < 1e-8, break, end
    cms.relax;
end
fprintf('  Density profile fixed.\n')

%% Compare computed and analytic density structure
q = cms.qrot;
s3 = cms.bs(1); % mean radius, s^3=b*a^2 but a=1
m = q*s3;

% Zharkov & Trubistyn (1978) Table 3.1
ZT5 = [nan; 0.0830; 1.4798; 5.929; 3.497; 2.52; 2.4; nan; nan];
% Hubbard (2013) Table 1
H13_128 = [q; 0.082999915; 1.4798138; 5.9269129; 3.4935680; 2.5493209;...
    2.1308951; 1.9564143; 1.9237724];

% CMSPlanet
CMP = [q; m; cms.Jn(2)*1e2; -cms.Jn(4)*1e4; cms.Jn(6)*1e5; -cms.Jn(8)*1e6;...
    cms.Jn(10)*1e7; -cms.Jn(12)*1e8; cms.Jn(14)*1e9];

% Make it a table
T = table(ZT5, H13_128, CMP);
T.Properties.RowNames = {'q','m','J2x10^2','-J4x10^4','J6x10^5','-J8x10^6',...
    'J10x10^7','-J12x10^8','J14x10^9'};

% Display
format long
format compact
fprintf('\n')
disp(T)
format

%% Save and deliver
% save('linear_density_model', 'cms', 'T')
% try
% sendmail('address','n=1 polytrope','','index1polytrope.mat')
% catch
% end
% !shutdown /s /t 30
% exit force
