function cmp = double_polytrope_w_core(N, x, zstrat)
%DOUBLE_POLYTROPE_W_CORE Model planet with two polytropes and a solid core.
%    DOUBLE_POLYTROPE_W_CORE(N, x) returns an N-layer CMSPlanet object with
%    three barotropes eos objects. The first is a polytrope defined by constant
%    x(1) and index x(2) and applied to layers 1:tind-1. The second is a
%    polytrope defined by constant x(3) and index x(4) and applied to layers
%    tind:N-1. Finally, layer N is assigned a barotropes.ConstDensity eos with
%    rho0=x(5). The transition between the polytropes is at layer tind, the
%    layer with (ai/a0) nearest to x(6). The core (layer N) is given a
%    normalized equatorial radius equal to x(7).
%
%    The default layer spacing is one of equal radius increments between a/a0=1
%    and a/a0=1/N, except that the transition layers are "snapped" to the given
%    values, x(6) and x(7).
%
%    DOUBLE_POLYTROPE_W_CORE(N, x, zstrat) lets you specify the layer spacing.
%    Pass a handle to a function that takes a single scalar integer (number of
%    layers) and returns a vector of that length with values in the interval
%    [x(7),1], for normalized layer radii. See +zvecs for examples.

if nargin == 0
    help('generators.double_polytrope_w_core')
    return
end
narginchk(2,3)
if ((nargin < 3) || isempty(zstrat)), zstrat = @(n)linspace(1, x(7), n); end
validateattributes(N,{'numeric'},{'positive','integer'},'','N',1)
validateattributes(x,{'numeric'},{'nonnegative','vector','numel',7},2)
validateattributes(zstrat,{'function_handle'},{},'','zstrat',3)
assert(x(6) > 0 && x(6) < 1, 'Envelope transition radius must be in (0,1).')
assert(x(7) > 0 && x(7) < 1, 'Core radius must be in (0,1).')
assert(x(7) < x(6), 'Core must be interior to lower envelope.')

% Generate and verify the layer distribution
zvec = zstrat(N);
assert(isnumeric(zvec) && isvector(zvec) && (numel(zvec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(zvec > 0) && all(zvec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(zvec(end-1) > x(7),...
    'Chosen layer spacing would put one or more layers below the core.')
if abs(zvec(end) - x(7)) > (2/N)
    warning('Check radial distribution: is the core where you wanted it?')
end

% Find and verify eos transition
[~, tind] = min(abs(zvec-x(6)));
if tind == 1
    warning('CMP:GEN',...
        'Transition too close to surface; first polytrope will have zero layers.')
end
if tind == N
    warning('CMP:GEN',...
        'Transition too close to core; second polytrope will have zero layers.')
end

% Snap-to-grid transition locations
zvec(tind) = x(6);
zvec(end) = x(7);

% Create the CMSPlanet object with assigned layer distribution
cmp = CMSPlanet;
cmp.ai = zvec;
cmp.rhoi = ones(N,1);

% Create and assign the eoss
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.ConstDensity(x(5));

cmp.eos = [repmat(eos1, tind - 1, 1);...
           repmat(eos2, N - tind, 1);...
           eos3];

% Initialize density at 1-bar values, just because
for k=1:N
    cmp.rhoi(k) = cmp.eos(k).density(1e5);
end

end
