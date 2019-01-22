function cmp = layered_planet(N, EOSs, rts, lamstrat, forcematch)
%LAYERED_PLANET Convenience generator for assigining EOSs to layers.
%    LAYERED_PLANET(N, EOSs, rts) returns an N-layer CMSPlanet object with
%    length(EOSs) eos objects. Transition from EOSs(i) to EOSs(i+1) is at layer
%    tind(i), the layer with (ai/a0) nearest rts(i). EOSs(1) is applied to layers
%    1:tind(1)-1, EOSs(2) is applied to layers tind(1):tind(2)-1, and so forth.
%    EOSs(end) is applied to layers tind(end):N. Thus length(rts) must equal
%    length(EOSs)-1. The default layer spacing is the one returned by
%    lambdas.best(N).
%
%    LAYERED_PLANET(N, EOSs, rts, lamstrat) lets you specify the lambda spacing.
%    Pass a handle to a function that takes a single scalar integer (number of
%    layers) and returns a vector of that length with values in the interval (0,
%    1], for normalized layer radii. For example, to set layers with equally
%    spaced radii use lamstrat=@(n)linspace(1,1/n,n). Note that the final layer
%    radii might be slightly different due to placement of transition radii. A
%    collection of pre-made distributions is available in package +lambdas.
%
%    LAYERED_PLANET(..., forcematch) if forcematch=true forces the normalized
%    radii of the transition from EOSs(i) to EOSs(i+1) to exactly match rts(i).
%    This is applied after the initial lambda spacing. The default is
%    forcematch=true.
%
% Examples:
%    % planet with two-polytrope envelopes and constant density core
%    eoss = [barotropes.Polytrope(1e5, 1),...
%            barotropes.Polytrope(8e4, 1),...
%            barotropes.ConstDensity(8000)];
%    cmp = generators.layered_planet(128, eoss, [0.8, 0.2]);

% Handle inputs
if nargin == 0
    help('generators.layered_planet')
    return
end
narginchk(3,5)
if ((nargin < 4) || isempty(lamstrat)), lamstrat = @lambdas.best; end
if ((nargin < 5) || isempty(forcematch)), forcematch = true; end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(EOSs, {'barotropes.Barotrope'}, {'vector'}, '', 'EOSs', 2)
validateattributes(rts, {'numeric'}, {'vector'}, '', 'rts', 3)
validateattributes(lamstrat, {'function_handle'}, {}, '', 'lamstrat', 4)
validateattributes(forcematch, {'logical'}, {'scalar'}, '', 'forcematch', 5)
assert(all(rts > 0 & rts < 1), 'Transition (normalized) radii must be in (0,1).')
assert(isequal(rts, sort(rts,'descend')), 'Specify transition radii from top.')
assert(length(rts) == length(EOSs) - 1, 'Specify n-1 radii for n EOSs.')

% Create desired lambdas spacing
lams = lamstrat(N);
assert(isnumeric(lams) && isvector(lams) && (numel(lams) == N),...
    '@lamstrat(N) must return a vector of length N with values in (0,1].')
assert(all(lams > 0) && all(lams <= 1),...
    '@lamstrat(N) must return a vector of length N with values in (0,1].')

% Find transition radii and indices
tind = nan(size(rts));
for k=1:length(rts)
    [~, tind(k)] = min(abs(lams-rts(k)));
end
assert(tind(1) > 1,...
    'First transition too close to surface; first eos has zero layers.')
for k=1:length(tind)-1
    assert(tind(k+1) > tind(k),...
        'Transitions are too close together; eos(%g) has zero layers.', k)
end
assert(tind(end) < N,...
    'Last transition too close to center; last eos has zero layers.')
if forcematch, lams(tind) = rts; end

% Construct new planet and assign lambdas and eoss
cmp = CMSPlanet;
cmp.ai = lams;
cmp.rhoi = ones(N,1);
cmp.eos = repmat(barotropes.ConstDensity(0), N, 1);
for k=1:N
    zind = find(rts < lams(k), 1);
    if ~isempty(zind)
        cmp.eos(k) = EOSs(zind);
    else
        cmp.eos(k) = EOSs(end);
    end;
end

% For show, set a 1-bar density throughout (helps visualize the eos layers)
for k=1:N
    cmp.rhoi(k) = cmp.eos(k).density(1e5);
end

end
