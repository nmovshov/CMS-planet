% +ZVECS CMS layer distributions.
%
% Each function in this package returns a vector of length N suitable for use
% as the zvec argument to cms.m or, more typically, to create the ai property
% of a CMSPlanet object. For example:
%    cmp = CMSPlanet();
%    cmp.radius = 71492*1e3;
%    cmp.si = cmp.radius*zvecs.topheavy(1024);
%
% Distributions:
%   best         - Return our current notion of the best zvec distribution.
%   cosine       - Return a zvec distribution with cosine-like spacing.
%   equal_dr     - Return a zvec distribution with equally spaced radii.
%   equal_volume - Return a zvec distribution making layers of equal volume.
%   exponential  - Return a zvec distribution with exponentially spaced radii.
%   topheavy     - Return a zvec distribution with top-heavy spacing.
%   toptopheavy  - Return a zvec distribution with top-top-heavy spacing.
%   trizone      - Return a 3-zone zvec distribution.
