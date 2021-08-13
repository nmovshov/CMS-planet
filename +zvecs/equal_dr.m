function zvec = equal_dr(N, halftop)
%EQUAL_DR Return a zvec distribution with equally spaced radii.
%    zvec = EQUAL_DR(N) returns an N-vector of equally spaced normalized radii.
%
%    zvec = EQUAL_DR(N, halftop) where halftop==true makes the thickness of the
%    first layer be exactly half that of the other layers. The default is
%    halftop=true. I still do not understand this mystery hack.

narginchk(1,2)
if nargin < 2, halftop = true; end
validateattributes(N,{'numeric'},{'positive','integer','scalar'})
validateattributes(halftop,{'logical'},{'scalar'})

dl = 1/N;
if halftop
    zvec = [1, (1 - dl/2):-dl:1/N];
else
    zvec = 1:-dl:1/N;
end

end
