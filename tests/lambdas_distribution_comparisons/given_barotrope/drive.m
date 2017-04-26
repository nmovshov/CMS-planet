%% test convergence 
clear
clc

f = @(n)lambdas.topheavy(n,[3/4,1/2],true); 
name = 'top heavy [3/4, 1/2]';

x = 2.^[7:10];
for k=1:length(x)
    cmp = models.triple_polytrope(x(k), [1.14e5, 1, 4.6e4, 0.9, 30, 0.67, 0.85, 0.1], f, true);
    cmp.M = 2e27;
    cmp.a0 = 7e7;
    cmp.qrot = 0.1;
    cmp.opts.verbosity = 2;
    cmp.opts.MaxIterBar = 60;
    
    cmp.relax_to_barotrope;
    L(k) = cmp.NMoI;
    beta(k) = cmp.betanorm;
    y2(k) = cmp.J2;
    y4(k) = cmp.J4;
    y6(k) = cmp.J6;
    clear cmp
end
clear k
