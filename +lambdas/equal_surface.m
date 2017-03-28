%% EQUISURFACE LAMBA DISTRIBUTION
% Return a lambda distribution whose spheroids have equal surface.
function lambdas = equi_surface(N)

lambdas = power((1:N)/N, 1/2);

lambdas = flip(lambdas);

end