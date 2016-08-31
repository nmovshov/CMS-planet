%% CMS CONVERGENCE TESTS
% This script runs some numerical tests designed to verify convergence and help
% us choose good parameters.

%% Prepare workspace
clear
clc
close all
addpath(fullfile(pwd,'..'))

%% Test convergence of Laplace expansion

% Construct a CMS object and with a linear density profile
cms = ConcentricMaclaurinSpheroids;
cms.nlayers = 12;
cms.qrot = 0.1;
cms.deltas = linspace(0,1000,cms.nlayers); %TODO: that's not a linear profile

% Relax to hydrostatic equilibrium and validate
cms.relax();
cms.validate();
