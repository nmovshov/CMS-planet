% +LAMBDAS CMS layer distributions.
%
% Each function in this package creates a vector of length N suitable for use as
% the lambdas property of the ConcentricMaclaurinSpheroids class. For example:
%    cmp = CMSPlanet(128);
%    cmp.cms.lambdas = lambdas.topheavy(128);
%
% Files
%   cosine   - Return a lambda distribution with cosine-like spacing.
%   equal_dr - Return a lambda distribution with equally spaced radii.
%   topheavy - Return a lambda distribution with top-heavy spacing.
%   best     - Return our current notion of the best lambda distribution.
%   trizone  - Return a 3-zone lambda distribution.

