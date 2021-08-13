function zvec = equal_volume(N, halftop)
%EQUAL_VOLUME Return a zvec distribution making layers of equal volume.
%    zvec = EQUAL_VOLUME(N) returns an N-vector of normalized radii such that
%    the volume of the layer between zvec(i) and zvec(i+1) is constant.

narginchk(1,2)
if nargin < 2, halftop = true; end
validateattributes(N,{'numeric'},{'positive','integer','scalar'})
validateattributes(halftop,{'logical'},{'scalar'})

zvec = ((N:-1:1)/N).^(1/3);
if halftop
    dl = 1 - zvec(2);
    zvec(2:end) = zvec(2:end) + dl/2;
end

end
