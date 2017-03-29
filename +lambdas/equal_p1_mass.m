function lambdas = equal_p1_mass(N, halftop)
%EQUAL_P1_MASS Return a lambda distribution with approximately equal mass layers.
%    lambdas = EQUAL_P1_MASS(N) returns an N-vector of normalized radii making
%    layers of equal mass if the density follows an index-1 polytropic relation.
%
%    lambdas = EQUAL_P1_MASS(N, halftop) where halftop==true modifies the 2nd
%    element of lambdas to make the thickness of the first layer be exactly half
%    that of the second layer. The default is halftop=false.

narginchk(1,2)
if nargin == 1, halftop = false; end
validateattributes(halftop, {'logical'}, {'scalar'})

mr = @(ri)(sin(pi*ri) - pi*ri.*cos(pi*ri))/pi;
lambdas = zeros(1, N);
for k = 1:N;
    lambdas(k) = fzero(@(ri)(mr(ri) - k/N), k/N);
end
lambdas = flip(lambdas);

if halftop
    lambdas(2) = 1 - (1 - lambdas(3))/3;
end

end
