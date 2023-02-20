%% How many density layers and shape layers are required for desired precision?
clear
clc
close all

%% Choose layer numbers and distribution strategies to investigate
N = [128, 256];
nx = [32, 64];

%% The experiment MAY be sensitive to the chosen density model
% Here we use a smooth polynomial model, with no discontinuities or even steep
% slopes. The goal is to find the *minimal* resolution needed for a desired J
% precisoion. More complex density profiles will surely require higher
% resolution, and/or special treatment of discontinuous density.

rho0 = 0.199;
p=[-9.1352277e+02, -2.0616194e+03, -2.5768977e+02, -6.8877550e+01,...
    8.6817818e+03, 4.2076235e+03, -1.3579913e+04, 0.0, 3.9924165e+03];
qrot = 0.1574;

%% Create a table to store experiment results
vars = {'N','nx','runtime','J2','J4'};
tps = {'double','double','double','double','double'};
nN = length(N);
nnx = length(nx);
T = table('Size', [nN*nnx, length(vars)],...
    'VariableTypes', tps,...
    'VariableNames', vars);

%% Call cms, this will take a while!!
fname = sprintf('%f.mat',now());
fprintf("Running CMS, this will take a while!\n")
for j=1:nN
    for k=1:nnx
        row = (j-1)*(nnx) + k;
        fprintf('working on row %d of %d (N=%d, nx=%d)...',row,height(T),N(j),nx(k));
        T.N(row) = N(j);
        T.nx(row) = nx(k);
        zvec = linspace(1, 0.01, N(j));
        dvec = polyval(p, zvec);
        cp = CMSPlanet.default_planet();
        cp.opts.xlayers = nx(k);
        cp.ai = zvec;
        cp.rhoi = dvec;
        cp.optimize_lambdas();
        tic
        cp.relax_to_HE();
        T.runtime(row) = toc;
        T.J2(row) = cp.J2;
        T.J4(row) = cp.J4;
        fprintf('done.\n')
        save(fname);
    end
end
fprintf('All done. (%s)\n',seconds2human(sum(T.runtime)));
