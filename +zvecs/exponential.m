function zvec = exponential(N, halftop)
%EXPONENTIAL Return a lambda distribution with exponentially spaced radii.
%    zvec = EXPONENTIAL(N) returns an N-vector of exponentially spaced
%    normalized radii.
%
%    zvec = EXPONENTIAL(N, halftop) where halftop==true makes the thickness of
%    the first layer be exactly half that of the other layers. The default is
%    halftop=true. I still do not understand this mystery hack.

narginchk(1,2)
if nargin < 2, halftop = true; end
validateattributes(N,{'numeric'},{'positive','integer','scalar'})
validateattributes(halftop,{'logical'},{'scalar'})

zvec = log((exp(2)-1)*(N:-1:1)/N+1)/2;

if nargin == 1, halftop = false; end
if halftop
    dl = 1 - zvec(2);
    zvec(2:end) = zvec(2:end) + dl/2;
end

end
