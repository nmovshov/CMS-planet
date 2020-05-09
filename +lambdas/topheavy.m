function lambdas = topheavy(N, skew, I, halftop)
%TOPHEAVY Return a lambda distribution with top-heavy spacing.
%    lambdas = TOPHEAVY(N) returns a vector of normalized radii with three
%    quarters of them equally distributed in the top half and the rest equally
%    distributed in the bottom half of the interval [1/N,1].
%
%    lambdas = TOPHEAVY(N, skew) returns a vector of normalized radii with skew(1)
%    of them equally distributed in the interval [1-skew(2), 1] and the rest
%    equally distributed in the interval [1/N, 1-skew(2)].
%
%    lambdas = TOPHEAVY(N, skew, interval) returns a vector of normalized radii
%    with skew(1) of them equally distributed in the interval
%    [I(1)+(I(2)-I(1))*(1-skew(2)),I(2)] and the rest equally distributed in the
%    interval (I(1),I(1)+(I(2)-I(1))*(1-skew(2))].
%
%    lambdas = TOPHEAVY(..., halftop) where halftop==true makes the thickness of
%    the first layer be exactly half that of the other layers in the upper zone.
%    The default is halftop=true.

if nargin == 0
    help('lambdas.topheavy')
    return
end
narginchk(1,4)
if nargin < 2 || isempty(skew), skew = [3/4, 1/2]; end
if nargin < 3 || isempty(I), I = [1/N,1]; end
if nargin < 4 || isempty(halftop), halftop = true; end
validateattributes(skew, {'numeric'}, {'vector', 'numel', 2, '>', 0, '<', 1})
validateattributes(I, {'numeric'}, {'nonnegative', 'vector', 'numel', 2})
validateattributes(halftop, {'logical'}, {'scalar'})

n1 = fix(skew(1)*N);
n2 = N - n1;
top = [I(1) + (I(2) - I(1))*(1 - skew(2)), I(2)];
bot = [I(1), I(1) + (I(2) - I(1))*(1 - skew(2))];
lam1 = linspace(top(2), top(1), n1);
lam2 = linspace(bot(2), bot(1), n2+1);
lambdas = flip(unique([lam1, lam2]));

if halftop
    dl = I(2) - lambdas(2);
    lambdas(2:end) = lambdas(2:end) + dl/2;
end

end
