%% How many density layers and shape layers are required for desired precision?
clear
clc
close all

%% Choose layer numbers and distribution strategies to investigate
N = [128, 256];
nx = [32, 64];
zstrat = {@zvecs.equal_dr};

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
vars = {'N','nx','zstrat','runtime','J2','J4'};
tps = {'double','double','string','double','double','double'};
nN = length(N);
nZ = length(zstrat);
nnx = length(nx);
T = table('Size', [nN*nnx*nZ, length(vars)],...
    'VariableTypes', tps,...
    'VariableNames', vars);

%% Call cms, this will take a while!!
fname = sprintf('%f.mat',now());
fprintf("Running CMS, this will take a while!\n")
for i=1:nZ
    for j=1:nN
        for k=1:nnx
            row = (i-1)*(nN*nnx) + (j-1)*(nnx) + k;
            fprintf('working on row %d of %d (N=%d, nx=%d)...',row,height(T),...
                N(j),nx(k));
            T.zstrat(row) = string(char(zstrat{i}));
            T.N(row) = N(j);
            T.nx(row) = nx(k);
            zvec = zstrat{i}(N(j));
            dvec = polyval(p, zvec);
            tic
            Js = cms(zvec, dvec, qrot, 'xlayers', nx(k));
            T.runtime(row) = toc;
            T.J2(row) = Js(2);
            T.J4(row) = Js(3);
            fprintf('done.\n')
            save(fname);
        end
    end
end
fprintf('All done. (%s)\n',seconds2human(sum(T.runtime)));
