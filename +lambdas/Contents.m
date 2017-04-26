% +LAMBDAS CMS layer distributions.
%
% Each function in this package returns a vector of length N suitable for use as
% the lambdas property of the ConcentricMaclaurinSpheroids class. For example:
%    cmp = CMSPlanet(128);
%    cmp.cms.lambdas = lambdas.topheavy(128);
%
% Distributions:
%   best         - Return our current notion of the best lambda distribution.
%   cosine       - Return a lambda distribution with cosine-like spacing.
%   equal_dr     - Return a lambda distribution with equally spaced radii.
%   equal_volume - Return a lambda distribution making layers of equal volume.
%   exponential  - Return a lambda distribution with exponentially spaced radii.
%   topheavy     - Return a lambda distribution with top-heavy spacing.
%   toptopheavy  - Return a lambda distribution with top-top-heavy spacing.
%   trizone      - Return a 3-zone lambda distribution.
