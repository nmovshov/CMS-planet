function cmp = single_polytrope(N, x)
%SINGLE_POLYTROPE The simplest toy model planet.
%    SINGLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object constructed with
%    a scalar barotropes.Polytrope eos using the polytropic constant x(1) and
%    the polytropic index x(2).

cmp = CMSPlanet(N);
cmp.eos = barotropes.Polytrope(x(1), x(2));

end
