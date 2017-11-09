%% test convergence of NMoI
clear
clc

f = @(n)lambdas.unused.trizone(n,[2/3,2/3],true); 
name = 'tri-zone [2/3,2/3] w/ halftop';

x = 2.^[7:11];
for k=1:length(x)
    cmp = cmsmodels.triple_polytrope(x(k), [1.14e5, 1, 4.6e4, 0.9, 30, 0.67, 0.85, 0.1], f, true);
    cmp.M = 2e27;
    cmp.a0 = 7e7;
    cmp.qrot = 0;
    cmp.opts.verbosity = 2;
    cmp.opts.MaxIterBar = 60;
    
    cmp.relax_to_barotrope;
    L(k) = cmp.NMoI;
    beta(k) = cmp.betanorm;
    clear cmp
end
clear k
