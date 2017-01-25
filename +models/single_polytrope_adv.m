function cmp = single_polytrope_adv(N, x)
%SINGLE_POLYTROPE_ADV The simplest toy model planet (advanced version).
%    SINGLE_POLYTROPE_ADV(N, x) returns an N-layer CMSPlanet object with a
%    barotropes.Polytrope with constant x(1) and index x(2) assigned to layers
%    2:N. Layer 1 is assigned a zero density barotropes.ConstDensity eos.

cmp = CMSPlanet(N);
eos0 = barotropes.ConstDensity(0);
eos1 = barotropes.Polytrope(x(1), x(2));
cmp.eos = [eos0; repmat(eos1, N-1, 1)];

end
