%% test convergence 
clear
clc

f = @(n)lambdas.topheavy(n,[3/4,2/5],true); 
name = 'top heavy [3/4,2/5] w/ halftop';

in = load('rho_of_r.mat');
x = 2.^[7:10];
for k=1:length(x)
    cmp = CMSPlanet(x(k));
    cmp.M = 2e27;
    cmp.a0 = 7e7;
    cmp.cms.lambdas = f(x(k));
    cmp.rhoi = in.rho_c*in.RHO(cmp.cms.lambdas);
    cmp.qrot = 0.1;
    cmp.opts.verbosity = 3;
    cmp.opts.MaxIterBar = 60;
    
    cmp.relax_to_HE;
    cmp.renormalize_density;
    L(k) = cmp.NMoI;
    beta(k) = cmp.betanorm;
    y2(k) = cmp.J2;
    y4(k) = cmp.J4;
    y6(k) = cmp.J6;
    clear cmp
end
clear k
