%% SPEEDTEST - ESTIMATE EXPECTED RUNTIME OF CMS MODELS

%% Prepare workspace
clear
clc
close all
addpath(fullfile(getenv('userprofile'),'CMS-planet'))

%% Time cms.relax() for several values of nlayers
% I will use a CMS model with a linear density profile and equally spaced
% layers. I give the planet fairly vigorous spin to ensure convergence will take
% a while. I measure the timing of 10 consecutive passes updating the shape
% functions and gravity moments (calculating shape functions dominates
% computation time) and save the mean time per pass.

nvec = [2, 4, 8, 16, 32, 64];
tvec = nan(size(nvec));

for k=1:length(nvec)
    cms = ConcentricMaclaurinSpheroids(nvec(k));
    cms.deltas = linspace(0,1,nvec(k));
    cms.lambdas = linspace(1,1/nvec(k),nvec(k)); % this is actually the default
    cms.qrot = 0.1;
    cms.opts.verbosity = 2;
    cms.opts.MaxIter = 10; % usually NOT enough to achieve convergence
    
    tvec(k) = cms.relax()/cms.opts.MaxIter;
end

%% Extrapolate timing based on asymptotic scaling of t~O(n^2)
p = polyfit(nvec.^2, tvec, 1);
N = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];
T = polyval(p,N.^2);

%% Display timing results and extrapolated predictions
ah = axes;
hold(ah, 'on')
plot(nvec, tvec, '-o', 'linewidth', 2, 'displayname', 'measured')
plot(N, T, '--', 'displayname', 'extrapolated')
ah.XScale = 'log';
ah.YScale = 'log';
xlabel('nlayers')
ylabel('time per pass [s]')
title(['CMS model timing with linear density profile and q = ',...
    num2str(cms.qrot)])
ah.Box = 'on';
legend show
