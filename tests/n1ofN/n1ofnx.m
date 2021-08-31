%% Test CMS convergence with increasing xlayers on the N=1 polytrope
% This script runs the n=1 polytrope with large N and varying xlayers to test
% convergence and runtime performance of skip-n-spline with CMSPlanet.
%
% Convergence is gauged against the full layer count case, but just for fun we
% aim for the numbers of Wisdom and Hubbard (2016) Table 3.

%% Prepare workspace
clear
clc
close all

%% Choose layer numbers and distribution strategies to investigate
N = 4096;
nx = [64, 128, 256, 512, 1024, 2048, 4096];
zstrat = @(n)zvecs.equal_dr(n);

%% Construct a polytrope of index 1, aiming for exact replicaiton of WH16
G = 6.6738480e-11; % Hubbard to Guillot to me personal communcation
GM = 1.266865361e17; % WH16
M = GM/G;
Re = 71492*1e3; % (to match K use K=2*G/pi*R^2 instead, but it doesn't matter)
qrot = 0.089195487; % WH16
wrot = sqrt(qrot*GM/Re^3);
Prot = 2*pi/wrot;
aos = 1.022875431133185; % WH16 Table 3 Re/R
K = 2.003565e5; % Hubbard to Guillot personal communication (no effect on Js)
n = 1;
eos = barotropes.Polytrope(K, n);
eos.name = '$P\propto\rho^2$';

%% Set up CMSPlanet(s)
for k=1:length(nx)
    cmp = CMSPlanet('xlayers',nx(k));
    cmp.name = sprintf('%d/%d-%s CMP',...
        nx(k),N,string(char(zstrat)).extractBetween('.','('));
    cmp.G = G; % undocumented TOFPlanet property
    cmp.mass = M;
    cmp.radius = Re;
    cmp.period = Prot; % trying to match WH16 qrot
    cmp.ai = Re*zstrat(N); % will be renormalized
    cmp.rhoi = ones(N,1)*M/(4*pi/3*Re^3); % will be renormalized
    cmp.P0 = 0.0; % added to surface pressure
    cmp.eos = eos;
    CMPS(k) = cmp;
end

%% Relax to desired barotrope (slow!)
fprintf('Running CMS (this will take a while!)...\n')
t = tic;
for k=1:length(nx)
    fprintf('Running case %d of %d (N=%d,nx=%d)...',k,length(nx),N,nx(k))
    cmp = CMPS(k);
    cmp.opts.verbosity = 0;
    cmp.opts.drhotol = 1e-6;
    cmp.opts.dJtol = 1e-8;
    cmp.opts.MaxIterBar = 60;
    cmp.opts.MaxIterHE = 2;
    rt(k) = cmp.relax_to_barotrope();
    fprintf(' done. (%s)\n', seconds2human(rt(k)))
end
fprintf('All done. (%s)\n',seconds2human(toc(t)))

%% Construct the benchmarking table
% The variables to compare are [Re/R, J2, J4, ..., J14]

% Wisdom and Hubbard (2016) Table 3
CLC = [nan, nan, nan, 1.022875431133185,...
        1.398851089834637e-2, -5.318281001092471e-4, 3.011832290533577e-5,...
       -2.132115710726158e-6, 1.740671195871128e-7, -1.568219505602588e-8,...
        1.518099230068580e-9];

% With CMSPlanet
for k=1:length(nx)
    cmp = CMPS(k);
    MCMS(k,:) = [N, nx(k), rt(k), cmp.a0/cmp.s0, cmp.Js(2:8)];
end

cols = {'N', 'nx', 'runtime', 'Re/R', 'J2', 'J4', 'J6', 'J8', 'J10', 'J12', 'J14'};
rows = {'CLC', CMPS.name};
A = [CLC; MCMS];
E = (A(2:end-1,:) - A(end,:))./A(end,:);
E(:,1:3) = A(2:end-1,1:3);
T_vals = array2table(A, 'VariableNames', cols, 'RowNames', rows);
T_errs = array2table(E, 'VariableNames', cols, 'RowNames', rows(2:end-1));

%% Output and save
format shorte
display(T_errs(:,4:7))
format
clear CMPS k n u
save(sprintf('%f.mat',now()), 'T_vals','T_errs')
