%% EQUISPACED LAMBA DISTRIBUTION
% Return a lambda distribution which is equispaced in radius.
function lambdas = equi_spaced(N)

lambdas = linspace(1/N, 1, N);

lambdas = flip(lambdas);

end