function zvec = toptopheavy(N, skew, halftop)
%TOPTOPHEAVY Return a zvec distribution with top-top-heavy spacing.
%    zvec = TOPTOPHEAVY(N) returns a vector of normalized radii with two thirds
%    of them distributed in the top half and the rest equally distributed in
%    the bottom half of the interval (0,1]. The upper distribution follow the
%    same rule: two thirds of those layers in the top half and the rest in the
%    bottom half of the interval [0.5,1].
%
%    zvec = TOPTOPHEAVY(N, skew) returns a vector of normalized radii with
%    skew(1) of them distributed in the interval [1-skew(2),1] and the rest
%    equally distributed in the interval (0,1-skew(2)]. The upper distribution
%    follow the same rule: skew(1) of those layers in the interval
%    [1-skew(2)^2,1] and the rest in the interval [1-skew(2),1-skew(2)^2].
%
%    zvec = TOPTOPHEAVY(..., halftop) where halftop==true makes the thickness
%    of the first layer be exactly half that of the other layers in the upper
%    zone. The default is halftop=true.

narginchk(1,3)
if (nargin < 2) || isempty(skew), skew = [2/3, 1/2]; end
if nargin < 3, halftop = true; end
validateattributes(N,{'numeric'},{'positive','integer','scalar'})
validateattributes(skew,{'numeric'},{'vector','numel',2,'>',0,'<',1})
validateattributes(halftop,{'logical'},{'scalar'})

n1 = fix(skew(1)^2*N);
n2 = fix(skew(1)*N-n1);
n3 = N - n1 - n2;
lam1 = linspace(1, 1 - skew(2)^2, n1);
lam2 = linspace(1 - skew(2)^2, 1 - skew(2), n2+1);
lam3 = linspace(1 - skew(2), 1/N, n3+1);
zvec = flip(unique([lam1, lam2, lam3]));

if halftop
    dl = 1 - zvec(2);
    zvec(2:end) = zvec(2:end) + dl/2;
end

end
