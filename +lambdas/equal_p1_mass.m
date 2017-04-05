function lambdas = equal_p1_mass(N)
%EQUAL_P1_MASS Return a lambda distribution with approximately equal mass layers.
%    lambdas = EQUAL_P1_MASS(N) returns an N-vector of normalized radii making
%    layers of equal mass if the density follows an index-1 polytropic relation.

mr = @(ri)(sin(pi*ri) - pi*ri.*cos(pi*ri))/pi;
lambdas = zeros(1, N);
for k = 1:N;
    lambdas(k) = fzero(@(ri)(mr(ri) - k/N), k/N);
end
lambdas = flip(lambdas);

end
