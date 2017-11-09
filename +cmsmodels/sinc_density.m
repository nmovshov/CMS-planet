function cmp = sinc_density(N, x, lamstrat)
%SINC_DENSITY Planet model with a density structure given by sinc.
%    SINC_DENSITY(N, x) returns an N-layer CMSPlanet object with a density
%    structure given by an expansion of sinc functions whose coefficients are in
%    x. The density will be:
%
%    rho(lamb) = x(1)*sinc(pi*lamb) + x(2)*sinc(2*pi*lamb) + ...
%
%    The degree of the expansion is selected to agree with the length of the
%    vector x given in input. Minimum length is one. The central density,
%    sum(x(i)), must be positive. The top layer will always have zero density.
%    Note that the actual rhoi quantity of output cmp model will be scaled to
%    match the mass and radius of the model.
%
%    SINC_DENSITY(N, x, lamstrat) where lamstrat is a 2-element vector lets you
%    specify the layer spacing strategy. Approximately lamstrat(1) of available
%    layers will be distributed in the top lamstrat(2) of the planet. For example,
%    passing [3/4, 0.2] concentrates the layers heavily in the top 20% of the
%    planet, leaving about N/4 layers to fill the bottom 80%. A single half-width
%    layer of zero density is still reserved for the surface. To use the default
%    spacing pass lambda=[].
%
%    SINC_DENSITY(N, x, lamstrat) where lamstrat is a function handle lets you
%    specify the lambda spacing completely. Pass a handle to a function that takes
%    a single scalar integer (number of layers) and returns a vector of that
%    length with values in the interval (0, 1], for normalized layer radii. For
%    example, to set layers with equally spaced radii use
%    lamstrat=@(n)linspace(1,1/n,n).

narginchk(2,3)
if ((nargin == 2) || isempty(lamstrat)), lamstrat = [3/4, 1/2]; end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(x, {'numeric'}, {'vector',}, '', 'x', 2)
assert(sum(x) > 0, 'The central density (sum(x(i))) must be positive.')
validateattributes(lamstrat, {'numeric','function_handle'}, {}, '', 'lamstrat', 3)
if isnumeric(lamstrat)
    validateattributes(lamstrat, {'numeric'},...
        {'vector', 'numel', 2, '>', 0, '<', 1}, '', 'lamstrat', 3)
end

cmp = CMSPlanet(N);

if (isa(lamstrat, 'function_handle'))
    lams = lamstrat(N);
    assert(isnumeric(lams) && isvector(lams) && (numel(lams) == N),...
        '@lamstrat(N) must return a vector of length N with values in (0,1].')
    assert(all(lams > 0) && all(lams <= 1),...
        '@lamstrat(N) must return a vector of length N with values in (0,1].')
    cmp.cms.lambdas = lams;
else
    n1 = fix(lamstrat(1)*(N - 1));
    n2 = N - n1 - 1;
    dl1 = lamstrat(2)/(n1 - 1);
    dl2 = (1 - lamstrat(2))/(n2 + 1);
    lam1 = linspace(1 - dl1/2, (1 - lamstrat(2)), n1);
    lam2 = linspace((1 - lamstrat(2)) - dl2, dl2, n2);
    cmp.cms.lambdas = [1, lam1, lam2]';
end

order = length(x);

rhoi = zeros(N, 1);
for n = 1:order;
    rhoi = rhoi + x(n)*sin(n*pi*cmp.cms.lambdas)./(n*pi*cmp.cms.lambdas);
end

cmp.cms.deltas = [rhoi(1); diff(rhoi)];

% Some warning signs
if any(rhoi < 0)
    warning('Some densities are negative!')
end
if any(cmp.cms.deltas(2:end) < 0)
    warning('The density structure is not monotonic!')
end

end
