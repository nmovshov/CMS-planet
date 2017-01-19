function cmp = single_polytrope_w_core(N, x)
%SINGLE_POLYTROPE_W_CORE Toy model planet approximated by a polytrope and core.
%    SINGLE_POLYTROPE_CORE(N, x) returns an N-layer CMSPlanet object constructed
%    with a barotropes.Polytrope eos using the polytropic constant x(1) and the
%    polytropic index x(2) applied to lower-indexed layers and a
%    barotropes.ConstDensity eos with density x(3) applied to higher-indexed
%    layers. Transition is at layer index fix(N*(1 - x(4))). Layer equatorial
%    radii are equally spaced.

cmp = CMSPlanet(N);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.ConstDensity(x(3));
tind = fix(N*(1 - x(4)));
cmp.eos = [repmat(eos1, tind, 1); repmat(eos2, N - tind, 1)];

end
