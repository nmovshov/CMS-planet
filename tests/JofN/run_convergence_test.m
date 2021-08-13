%% How many density layers and shape layers are required for desired precision?
clear
clc
close all

%% Choose layer numbers and lambda strategies to investigate
N = [128];
nx = [32];
lamstrat = {@lambdas.best};

%% The experiment MAY be mildly sensitive to the chosen density model
rho0 = 0.199;
p=[-9.1352277e+02 -2.0616194e+03 -2.5768977e+02 -6.8877550e+01 8.6817818e+03 4.2076235e+03 -1.3579913e+04 0.0 3.9924165e+03];
mdl = @(n,s)generators.polynomial(n,p,s);
qrot = 0.1574;

%% Create a table to store experiment results
vars = {'N','nx','lamstrat','J2','runtime','cmp'};
tps = {'double','double','string','double','double','CMSPlanet'};
T = table('Size', [length(N)*length(nx)*length(lamstrat), length(vars)],...
    'VariableTypes', tps,...
    'VariableNames', vars);
[~,fname,~] = fileparts(tempname);

%% Relax CMPs, this may take a while
fprintf("Relaxing CMPs, this may take a while!\n")
for i=1:length(lamstrat)
    for j=1:length(N)
        for k=1:length(nx)
            row = (i-1)*length(N) + (j-1)*length(nx) + k;
            fprintf('working on row %d of %d (N=%d, nx=%d)...',row,height(T),...
                N(j),nx(k));
            T.lamstrat(row) = string(char(lamstrat{i}));
            T.N(row) = N(j);
            T.nx(row) = nx(k);
            cmp = mdl(N(j), lamstrat{i});
            cmp.opts.xlayers = nx(k);
            cmp.qrot = qrot;
           trun = cmp.relax_to_HE;
            T.runtime(row) = trun;
            T.J2(row) = cmp.J2;
            T.cmp(row) = cmp;
            fprintf('done.\n')
            save;
        end
    end
end
fprintf('All done. (%s)\n',seconds2human(sum(T.runtime)));

%% Examine the convergence behavior
