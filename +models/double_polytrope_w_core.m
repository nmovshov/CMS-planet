function cmp = double_polytrope_w_core(N, x, lamstrat, forcematch)
%DOUBLE_POLYTROPE_W_CORE Two polytrope and core model planet.
%    DOUBLE_POLYTROPE_W_CORE(N, x) returns an N-layer CMSPlanet object with two
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to layers 2:tind-1. Second polytrope defined by
%    constant x(3) and index x(4) is applied to layers tind:N-1. Transition is at
%    layer tind, the layer with (ai/a0) nearest to x(5). Layer N is assigned a
%    barotropes.ConstDensity eos wirh rho0=x(6) and given a normalized equatorial
%    radius equal to x(7). Layer 1 is assigned a zero density
%    barotropes.ConstDensity eos and is approximately half the width of the layers
%    below. The default layer spacing concentrates 3/4 of the available layers in
%    the top 0.5 of the envelope.
%
%    DOUBLE_POLYTROPE_W_CORE(N, x, lamstrat) where lamstrat is a 2-element vector
%    lets you specify the layer spacing strategy. Approximately lamstrat(1) of
%    available layers will be distributed in the top lamstrat(2) of the envelope.
%    For example, passing [3/4, 0.2] concentrates the layers heavily in the top
%    20% of the envelope, leaving about N/4 layers to fill the bottom 80%. A
%    single half-width layer of zero density is still reserved for the surface. To
%    use the default spacing pass lambda=[].
%
%    DOUBLE_POLYTROPE_W_CORE(N, x, lamstrat) where lamstrat is a function handle
%    lets you specify the lambda spacing completely. Pass a handle to a function
%    that takes a single scalar integer (number of layers) and returns a vector of
%    that length with values in the interval [x(7), 1], for normalized layer
%    radii. For example, to set layers with equally spaced radii in the envelope
%    use lamstrat=@(n)linspace(1,x(7),n).
%
%    DOUBLE_POLYTROPE_W_CORE(..., forcematch) if forcematch=true forces the
%    normalized radius of the transition layer to exactly match x(5). This is
%    applied after the lambda spacing strategy. The default is forcematch=true.

narginchk(2,4)
if ((nargin == 2) || isempty(lamstrat)), lamstrat = [3/4, 1/2]; end
if ((nargin < 4) || isempty(forcematch)), forcematch = true; end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(x, {'numeric'}, {'vector', 'numel', 7, 'nonnegative'}, 2)
validateattributes(lamstrat, {'numeric','function_handle'}, {}, '', 'lamstrat', 3)
if isnumeric(lamstrat)
    validateattributes(lamstrat, {'numeric'},...
        {'vector', 'numel', 2, '>', 0, '<', 1}, '', 'lamstrat', 3)
end
validateattributes(forcematch, {'logical'}, {'scalar'}, '', 'forcematch', 4)
assert(x(5) > 0 && x(5) < 1, 'Transition (normalized) radius must be in (0,1).')
assert(x(7) > 0 && x(7) < 1, 'The core must have a normalized radius in (0,1).')
assert(x(7) < x(5), 'Core must be below transition radius.')

cmp = CMSPlanet(N);

if (isa(lamstrat, 'function_handle'))
    lams = lamstrat(N);
    assert(isnumeric(lams) && isvector(lams) && (numel(lams) == N),...
        '@lamstrat(N) must return a vector of length N with values in (0,1].')
    assert(all(lams > 0) && all(lams <= 1),...
        '@lamstrat(N) must return a vector of length N with values in (0,1].')
    lams = sort(lams, 'descend');
    assert(lams(end-1) > x(7),...
        'Chosen layer spacing would put one or more layers below the core.')
    if abs(lams(end) - x(7)) > (2/N)
        beep
        warning('Check lambda distribution: is the core where you wanted it?')
    end
    lams(end) = x(7);
    cmp.cms.lambdas = lams;
else
    n1 = fix(lamstrat(1)*(N - 1));
    n2 = N - n1 - 1;
    dl1 = lamstrat(2)*(1 - x(7))/(n1 - 1);
    dl2 = (1 - lamstrat(2))*(1 - x(7))/n2;
    lam1 = linspace(1 - dl1/2, 1 - lamstrat(2)*(1 - x(7)), n1);
    lam2 = linspace(1 - lamstrat(2)*(1 - x(7)) - dl2, x(7), n2);
    lams = [1, lam1, lam2]';
end

[~, tind] = min(abs(lams-x(5)));
assert(tind > 2,...
    'Transition too close to surface; first polytrope has zero layers.')
assert(tind < N,...
    'Transition too close to core; second polytrope has zero layers.')
if forcematch, lams(tind) = x(5); end

cmp.cms.lambdas = lams;

eos0 = barotropes.ConstDensity(0);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.ConstDensity(x(6));

cmp.eos = [eos0; repmat(eos1, tind - 2, 1); repmat(eos2, N - tind, 1); eos3];

end
