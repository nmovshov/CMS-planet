function lambdas = topheavy(N, skew, halftop)
%TOPHEAVY Return a lambda distribution with top-heavy spacing.
%    lambdas = TOPHEAVY(N) returns a vector of normalized radii with two thirds of
%    them equally distributed in the top half and the rest equally distributed in
%    the bottom half of the interval (0,1].
%
%    lambdas = TOPHEAVY(N, skew) returns a vector of normalized radii with skew(1)
%    of them equally distributed in the interval [1-skew(2), 1] and the rest
%    equally distributed in the interval (0, 1-skew(2)].
%
%    lambdas = TOPHEAVY(..., halftop) where halftop==true makes the thickness of
%    the first layer be exactly half that of the other layers in the upper zone.
%    The default is halftop=false.

narginchk(1,3)
if (nargin < 2) || isempty(skew), skew = [2/3, 1/2]; end
if nargin < 3, halftop = false; end
validateattributes(skew, {'numeric'}, {'vector', 'numel', 2, '>', 0, '<', 1})
validateattributes(halftop, {'logical'}, {'scalar'})

n1 = fix(skew(1)*N);
n2 = N - n1;
lam1 = linspace(1, 1 - skew(2), n1);
lam2 = linspace(1 - skew(2), 1/N, n2+1);
lambdas = flip(unique([lam1, lam2]));

if halftop
    dl = 1 - lambdas(2);
    lambdas(2:end) = lambdas(2:end) + dl/2;
end

end
