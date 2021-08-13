%% GRAVITY COEFFICIENTS OF ROTATING INDEX-1 POLYTROPE
% Example and test of the CMSPlanet class. We construct and converge a model of
% a rotating fluid planet with the pressure-density law
%
% $$P = K\rho^2$$
%
% with a polytropic constant $K$. This script demonstrates how to set up a
% CMSPlanet object with a specified barotrope and converge to a density
% structure in hydrostatic equilibrium with the given barotrope. The default
% starting density profile is that of a homogeneous sphere.
%
% The main output is a benchmarking table. The numbers to aim for are
% apparently those of the CLC method (Wisdom, 1996, unpublished) which is lost
% to time. I take the numbers from Wisdom and Hubbard (2016) Table 3. For
% consistency I try to replicate the Wisdom and Hubbard model exactly,
% including using their values for G and K (Guillot, double hearsay) and trying
% to match their q with a rotation period.
%
% Also comparing with CMS results from Wisdom and Hubbard (2016) Table 4, and
% with some q powers I dug out of the ZT red book that I think are supposed to
% match a 3rd-order ToF solution. Also adding some more numbers I got from
% Nadine (Nettelmann 2021 private communication).

%% Prepare workspace
clear
clc
close all

%% Construct a polytrope of index 1, aiming for exact replicaiton of WH16
G = 6.6738480e-11; % Hubbard to Guillot personal communcation
GM = 1.266865361e17; % WH16
M = GM/G;
Re = 71492*1e3; % (to match K use K=2*G/pi*R^2 instead)
qrot = 0.089195487; % WH16
wrot = sqrt(qrot*GM/Re^3);
Prot = 2*pi/wrot;
aos = 1.022875431133185; % WH16 Table 3 Re/R
K = 2.003565e5; % Hubbard to Guillot personal communication (no effect on Js)
n = 1;
eos = barotropes.Polytrope(K, n);
eos.name = '$P\propto\rho^2$';

%% Set up CMSPlanet
N = 128;
nx = 64;

cmp = CMSPlanet();
cmp.name = [int2str(N),'-point CMS'];
cmp.G = G; % undocumented CMSPlanet property
cmp.mass = M;
cmp.radius = Re;
cmp.period = Prot; % trying to match WH16 qrot
cmp.ai = Re*linspace(1, 1/N, N)'; % will be renormalized
cmp.rhoi = ones(N,1)*M/(4*pi/3*Re^3); % will be renormalized
cmp.P0 = 0; % added to surface pressure
cmp.eos = eos;

%% Relax to desired barotrope
cmp.opts.drhotol = 1e-6;
cmp.opts.dJtol = 1e-10;
cmp.opts.MaxIterBar = 60;
cmp.opts.MaxIterHE = 60;
cmp.opts.xlayers = nx;
cmp.relax_to_barotrope();

%% Construct the benchmarking table
% The variables to compare are [Re/R, J2, J4, ..., J14]
% With CMSPlanet
CMP = [cmp.a0/cmp.s0, cmp.Js(2:8)];
% Wisdom and Hubbard (2016) Table 3
CLC = [1.022875431133185, 1.398851089834637e-2, -5.318281001092471e-4,...
                          3.011832290533577e-5, -2.132115710726158e-6,...
                          1.740671195871128e-7, -1.568219505602588e-8,...
                          1.518099230068580e-9];
% Wisdom and Hubbard (2016) Table 4
CMS512 = [nan,            1.398924011471443e-2, -5.318792055591143e-4,...
                          3.012230402792236e-5, -2.132458660888379e-6,...
                          1.740988882124251e-7, -1.568529225158399e-8,...
                          1.518412099478797e-9];
% Hubbard (2013) Table 5
H13 =    [nan, 1.3989253e-2, -5.3187997e-4, 3.0122356e-5, -2.1324628e-6,...
               1.7409925e-7, -1.5685327e-8, 1.5184156e-9];

% Zharkov & Trubistyn (1978) eq. 34.12
q = cmp.qrot;
ZT3 = [nan, (0.173273*q - 0.197027*q^2 + 0.15*q^3),...
            (-0.081092*q^2 + 0.15*q^3),...
            (0.056329*q^3), nan, nan, nan, nan];

cols = {'Re/R', 'J2', 'J4', 'J6', 'J8', 'J10', 'J12', 'J14'};
rows = {'CLC', 'CMP', 'CMS_WH16', 'CMS_H13', 'ZT78'};
A = [CLC; CMP; CMS512; H13; ZT3];
E = (A(2:end,:) - A(1,:))./A(1,:);
T_vals = array2table(A, 'VariableNames', cols, 'RowNames', rows);
T_errs = array2table(E, 'VariableNames', cols, 'RowNames', rows(2:end));

%% Output
format shorte
%display(T_vals)
display(T_errs(:,1:3))
format
