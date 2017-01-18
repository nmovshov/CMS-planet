function cmp = double_polytrope(N, x)
%DOUBLE_POLYTROPE Toy model planet approximated by two polytropes.
%    DOUBLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object constructed with
%    two barotropes.Polytrope eos objects. First polytrope defined by constant
%    x(1) and index x(2), second polytrope defined by constant x(3) and index
%    x(4). Transition at layer index fix(N*(1 - x(5))).

cmp = CMSPlanet(N);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
tind = fix(N*(1 - x(5)));
cmp.eos = [repmat(eos1, tind, 1); repmat(eos2, N - tind, 1)];

end
