%% CMS CONVERGENCE TESTS
% This script runs some numerical tests designed to verify convergence and help
% us choose good parameters.

%% Prepare workspace
clear
clc
close all
addpath(fullfile(pwd,'..'))

%% Test convergence of Laplace expansion

% Construct a CMS object and give it a linear density profile
cms = ConcentricMaclaurinSpheroids;
cms.opts.nlayers = 12;
