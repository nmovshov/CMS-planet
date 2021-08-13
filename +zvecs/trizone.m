function zvec = trizone(N, parts, halftop)
%TRIZONE Return a 3-zone zvec distribution.
%    zvec = TRIZONE(N) returns a vector of normalized radii with two thirds of
%    them equally distributed in the top third, two thirds of the remainder
%    equally distributed in the middle third, and the rest equally distributed
%    in the bottom third of the interval (0,1].
%
%    zvec = TRIZONE(N, parts) returns a vector of normalized radii with
%    parts(1) of them equally distributed in the top third, parts(2) of the
%    reaminder equally distributed in the middle third, and the rest equally
%    distributed in the bottom third of the interval (0,1]. For example, the
%    default spacing above is obtained with parts=[2/3, 2/3].
%
%    zvec = TRIZONE(..., halftop) where halftop==true makes the thickness of
%    the first layer be exactly half that of the other layers in the upper
%    zone. The default is halftop=true.

narginchk(1,3)
if (nargin < 2) || isempty(parts), parts = [2/3, 2/3]; end
if nargin < 3, halftop = true; end
validateattributes(N,{'numeric'},{'positive','integer','scalar'})
validateattributes(parts,{'numeric'},{'vector','numel',2,'>=',0,'<',1})
validateattributes(halftop,{'logical'},{'scalar'})

n1 = fix(parts(1)*N);
n2 = fix(parts(2)*(N - n1));
n3 = N - n2 - n1;
lam1 = linspace(1, 2/3, n1);
lam2 = linspace(2/3, 1/3, n2+1);
lam3 = linspace(1/3, 1/N, n3+1);
zvec = flip(unique([lam1, lam2, lam3]));

if halftop
    dl = 1 - zvec(2);
    zvec(2:end) = zvec(2:end) + dl/2;
end

end
