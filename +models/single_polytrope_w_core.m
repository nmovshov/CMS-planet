function cmp = single_polytrope_w_core(N, x, lamstrat)
%SINGLE_POLYTROPE_W_CORE Polytrope-and-core planet.
%    SINGLE_POLYTROPE_W_CORE(N, x) returns an N-layer CMSPlanet object with a
%    barotropes.Polytrope with constant x(1) and index x(2) assigned to layers
%    2:N-1 and a barotropes.ConstDensity eos with rho0=x(3) assigned to layer N
%    which has a normalized equatorial radius equal to x(4). Layer 1 is assigned a
%    zero density barotropes.ConstDensity eos and is approximately half the width
%    of the layers below. The default layer spacing concentrates 2/3 of the
%    available layers in the top 0.5 of the envelope.
%
%    SINGLE_POLYTROPE_W_CORE(N, x, lamstrat) where lamstrat is a 2-element vector
%    lets you specify the layer spacing strategy. Approximately lamstrat(1) of
%    available layers will be distributed in the top lamstrat(2) of the envelope.
%    For example, passing [3/4, 0.2] concentrates the layers heavily in the top
%    20% of the envelope, leaving about N/4 layers to fill the bottom 80%. A
%    single half-width layer of zero density is still reserved for the surface. To
%    use the default spacing pass lambda=[].
%
%    SINGLE_POLYTROPE_W_CORE(N, x, lamstrat) where lamstrat is a function handle
%    lets you specify the lambda spacing completely. Pass a handle to a function
%    that takes a single scalar integer (number of layers) and returns a vector of
%    that length with values in the interval [x(4), 1], for normalized layer
%    radii. For example, to set layers with equally spaced radii in the envelope
%    use lamstrat=@(n)linspace(1,x(4),n).

narginchk(2,3)
if ((nargin == 2) || isempty(lamstrat)), lamstrat = [2/3, 1/2]; end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(x, {'numeric'}, {'vector', 'numel', 4, 'nonnegative'}, 2)
validateattributes(lamstrat, {'numeric','function_handle'}, {}, '', 'lamstrat', 3)
if isnumeric(lamstrat)
    validateattributes(lamstrat, {'numeric'},...
        {'vector', 'numel', 2, '>', 0, '<', 1}, '', 'lamstrat', 3)
end
assert(x(4) > 0 && x(4) < 1, 'The core must have a normalized radius in (0,1).')

cmp = CMSPlanet(N);

if (isa(lamstrat, 'function_handle'))
    lambdas = lamstrat(N);
    assert(isnumeric(lambdas) && isvector(lambdas) && (numel(lambdas) == N),...
        '@lamstrat(N) must return a vector of length N with values in (0,1].')
    assert(all(lambdas > 0) && all(lambdas <= 1),...
        '@lamstrat(N) must return a vector of length N with values in (0,1].')
    lambdas = sort(lambdas, 'descend');
    assert(lambdas(end-1) > x(4),...
        'Chosen layer spacing would put one or more layers below the core.')
    if abs(lambdas(end) - x(4)) > (2/N)
        beep
        warning('Check lambda distribution: is the core where you wanted it?')
    end
    lambdas(end) = x(4);
    cmp.cms.lambdas = lambdas;
else
    n1 = fix(lamstrat(1)*(N - 1));
    n2 = N - n1 - 1;
    dl1 = lamstrat(2)*(1 - x(4))/(n1 - 1);
    dl2 = (1 - lamstrat(2))*(1 - x(4))/n2;
    lam1 = linspace(1 - dl1/2, 1 - lamstrat(2)*(1 - x(4)), n1);
    lam2 = linspace(1 - lamstrat(2)*(1 - x(4)) - dl2, x(4), n2);
    cmp.cms.lambdas = [1, lam1, lam2]';
end

eos0 = barotropes.ConstDensity(0);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.ConstDensity(x(3));
cmp.eos = [eos0; repmat(eos1, N - 2, 1); eos2];

end
