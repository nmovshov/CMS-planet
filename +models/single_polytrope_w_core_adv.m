function cmp = single_polytrope_w_core_adv(N, x)
%SINGLE_POLYTROPE_W_CORE_ADV Model Polytrope-and-core model planet (advanced).
%    SINGLE_POLYTROPE_W_CORE_ADV(N, x) returns an N-layer CMSPlanet object 
%    with a barotropes.Polytrope with constant x(1) and polytropic index x(2)
%    assigned to layers 2:N-1 and a barotropes.ConstDensity eos with rho0 = x(3)
%    assigned to layer N, which has a normalized equatorial radius equal to x(4).
%    Layer 1 is assigned a zero density barotropes.ConstDensity eos. The layer
%    spacing is designed to minimize discretization error by concentrating 2/3 of
%    the available layers in the top 0.5 of the planet (see Hubbard & Militzer, 
%    2016).

cmp = CMSPlanet(N);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.ConstDensity(x(3));
cmp.cms.lambdas = linspace(1, x(4), N);
cmp.eos = [repmat(eos1, N-1, 1); eos2];

end
