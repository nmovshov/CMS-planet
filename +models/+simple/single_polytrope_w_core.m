function cmp = single_polytrope_w_core(N, x)
%SINGLE_POLYTROPE_W_CORE Toy model planet approximated by a polytrope and core.
%    SINGLE_POLYTROPE_W_CORE(N, x) returns an N-layer CMSPlanet object
%    constructed with a barotropes.Polytrope eos using the polytropic constant
%    x(1) and the polytropic index x(2) applied to the first N-1 layers and a
%    barotropes.ConstDensity eos with density x(3) applied to last layer. Layer
%    equatorial radii are equally spaced with the normalized radius of the last
%    layer (the "core") set to x(4).

cmp = CMSPlanet(N);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.ConstDensity(x(3));
cmp.cms.lambdas = linspace(1, x(4), N);
cmp.eos = [repmat(eos1, N-1, 1); eos2];

end
