%% Test CMS convergence with N on the n=1 polytrope
% This script benchmarks CMS with CMSPlanet against the CLC solution of WH16.
% The n=1 polytrope is run with cms with increasing number of layers. I use
% skip-n-spline conservatively.
%
% The numbers to benchmark against are apparently those of the CLC method
% (Wisdom, 1996, unpublished) which is lost to time. I take the numbers from
% Wisdom and Hubbard (2016) Table 3. For consistency we try to replicate the
% Wisdom and Hubbard model exactly, matching their q value exactly by using
% also their G and K values (Guillot, double hearsay).

%% Prepare workspace
clear
clc
close all

%% Choose layer numbers and distribution strategies to investigate
N = 2.^(8:13);
nx = -1;
zstrat = @(n)zvecs.topheavy(n);

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

%% Set up the CMSPlanet(s)
for k=1:length(N)
    cmp = CMSPlanet();
    cmp.name = sprintf('%d-%s CMP',...
        N(k),string(char(zstrat)).extractBetween('.','('));
    cmp.G = G; % undocumented CMSPlanet property
    cmp.mass = M;
    cmp.radius = Re;
    cmp.period = Prot; % trying to match WH16 qrot
    cmp.ai = Re*zstrat(N(k)); % will be renormalized
    cmp.rhoi = ones(N(k),1)*M/(4*pi/3*Re^3); % will be renormalized
    cmp.P0 = 0.0; % added to surface pressure
    cmp.eos = eos;
    CMPS(k) = cmp;
end

%% Relax to desired barotrope (slow!)
fprintf('Running CMS (this will take a while!)...\n')
t = tic;
for k=1:length(N)
    fprintf('Running case %d of %d (N=%d,nx=%d)...',k,length(N),N(k),nx)
    cmp = CMPS(k);
    cmp.opts.xlayers = nx;
    cmp.opts.verbosity = 0;
    cmp.opts.drhotol = 1e-6;
    cmp.opts.dJtol = 1e-8;
    cmp.opts.MaxIterBar = 60;
    cmp.opts.MaxIterHE = 2;
    rt(k) = cmp.relax_to_barotrope();
    fprintf(' done. (%s)\n', seconds2human(rt(k)))
    save
end
fprintf('All done. (%s)\n',seconds2human(toc(t)))

%% Construct the benchmarking table
% The variables to compare are [Re/R, J2, J4, ..., J14]

% Wisdom and Hubbard (2016) Table 3
CLC = [nan, nan, 1.022875431133185,...
        1.398851089834637e-2, -5.318281001092471e-4, 3.011832290533577e-5,...
       -2.132115710726158e-6, 1.740671195871128e-7, -1.568219505602588e-8,...
        1.518099230068580e-9];

% With CMSPlanet
for k=1:length(N)
    cmp = CMPS(k);
    MCMS(k,:) = [N(k),nx,cmp.a0/cmp.s0,cmp.Js(2:8)];
end

cols = {'N', 'nx', 'Re/R', 'J2', 'J4', 'J6', 'J8', 'J10', 'J12', 'J14'};
rows = {'CLC', CMPS.name};
A = [CLC; MCMS];
E = (A(2:end,:) - A(1,:))./A(1,:);
E(:,1:2) = A(2:end,1:2);
T_vals = array2table(A, 'VariableNames', cols, 'RowNames', rows);
T_errs = array2table(E, 'VariableNames', cols, 'RowNames', rows(2:end));

%% Output and save
format shorte
display(T_errs(:,3:4))
format
clear CMPS k n u
save(sprintf('%f.mat',now()), 'T_vals','T_errs')
