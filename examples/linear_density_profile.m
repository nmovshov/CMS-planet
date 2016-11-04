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
addpath(fullfile(pwd,'..'))
try
    si = setUnits; % if you have physunits
catch
    si = setFUnits; % if you don't have physunits
end
G = si.gravity;

%% Set up a CMS object
cms = ConcentricMaclaurinSpheroids(16);
cms.qrot = 0.088822426; % Hubbard (2013) Table 1

%% Construct a density profile linear in a_i to start
% A constant lambda step (default) calls for constant delta step
cms.deltas = ones(cms.nlayers, 1);
cms.deltas(1) = 0;

%% Relax to hydrostatic equilibrium
cms.opts.verbosity = 2;
cms.opts.MaxIterHE = 12;
cms.opts.email = '';
cms.relax;

%% After initial relaxation, we iteratively converge deltas to linear in s
iter = 0;
while (iter < 10)
    new_deltas = [0; -diff(cms.ss/cms.ss(1))];
    new_deltas = new_deltas/max(new_deltas); % just for show
    delta_deltas = max(abs(new_deltas - cms.deltas));
    fprintf('\nAdjusting density profile by < %g.\n',delta_deltas)
    if delta_deltas < 1e-2, break, end
    cms.deltas = new_deltas;
    iter = iter + 1;
    cms.relax;
end

%% Compare computed and analytic density structure
% q = cms.qrot;
% % Zharkov & Trubistyn (1978) eq. 34.12
% ZT3 = [q;...
%     (0.173273*q - 0.197027*q^2 + 0.15*q^3)*1e2;...
%     (-0.081092*q^2 + 0.15*q^3)*-1e4;...
%     (0.056329*q^3)*1e5;...
%     nan; nan; nan];
% % Hubbard (2013) Table 5
% H13_256 = [q; 1.3991574; 5.3182810; 3.0118323; 2.1321157; 1.7406710; 1.5682179];
% H13_512 = [q; 1.3989253; 5.3187997; 3.0122356; 2.1324628; 1.7409925; 1.5685327];
% 
% % CMSPlanet
% CMP = [q; cms.Js(2)*1e2; -cms.Js(3)*1e4; cms.Js(4)*1e5; -cms.Js(5)*1e6;...
%     cms.Js(6)*1e7; -cms.Js(7)*1e8;];
% 
% % Make it a table
% T = table(ZT3, H13_256, H13_512, CMP);
% T.Properties.RowNames = {'q','J2x10^2','-J4x10^4','J6x10^5','-J8x10^6',...
%     'J10x10^7','-J12x10^8'};
% 
% % Display
% format long
% format compact
% disp(T)
% format

%% Save and deliver
% save('index1polytrope', 'cmp', 'T')
% try
% sendmail('address','n=1 polytrope','','index1polytrope.mat')
% catch
% end
% !shutdown /s /t 30
% exit force
