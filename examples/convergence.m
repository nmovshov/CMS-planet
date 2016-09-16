%% CMS CONVERGENCE TESTS
% This script runs some numerical tests designed to verify convergence and help
% us choose good parameters.

%% Prepare workspace
clear
clc
close all
addpath(fullfile(pwd,'..'))

%% Test convergence of Laplace expansion

% Construct a CMS object with a linear density profile
N = 8;
cms = ConcentricMaclaurinSpheroids(N);
cms.qrot = 0.1;
cms.deltas = [0, ones(1,N-1)];

% Relax to hydrostatic equilibrium and validate
cms.relax();
cms.validate();
