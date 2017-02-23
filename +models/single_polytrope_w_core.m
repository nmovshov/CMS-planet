function cmp = single_polytrope_w_core(N, x)
%SINGLE_POLYTROPE_W_CORE Polytrope-and-core planet.
%    SINGLE_POLYTROPE_W_CORE(N, x) returns an N-layer CMSPlanet object with
%    a barotropes.Polytrope with constant x(1) and index x(2) assigned to layers
%    2:N-1 and a barotropes.ConstDensity eos with rho0=x(3) assigned to layer N
%    which has a normalized equatorial radius equal to x(4). Layer 1 is assigned
%    a zero density barotropes.ConstDensity eos and is approximately half the
%    width of the layers below. The layer spacing is designed to minimize
%    discretization error by concentrating 2/3 of the available layers in the
%    top 0.5 of the planet (see Hubbard & Militzer, 2016). However if x(4)>0.2
%    then layers are equally spaced above x(4).

cmp = CMSPlanet(N);

if x(4) > 0.2
    dl = (1 - x(4))/N;
    cmp.cms.lambdas = [1, linspace(1 - dl/2, x(4), N - 1)]';
else
    if x(4) <= 0, x(4) = eps; end
    n1 = fix(2/3*N) - 1;
    n2 = N - n1 - 1;
    dl1 = 0.5/(n1 - 1);
    dl2 = (0.5 - x(4))/n2;
    lam1 = linspace(1 - dl1/2, 0.5, n1);
    lam2 = linspace(0.5 - dl2, x(4), n2);
    cmp.cms.lambdas = [1, lam1, lam2]';
end

eos0 = barotropes.ConstDensity(0);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.ConstDensity(x(3));
cmp.eos = [eos0; repmat(eos1, N - 2, 1); eos2];

end
