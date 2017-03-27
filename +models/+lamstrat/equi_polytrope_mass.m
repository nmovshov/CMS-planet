%% EQUIMASS LAMBA DISTRIBUTION FOR POLYTROPE N=1
% Return a lambda distribution whose spheroids have equal mass for a
% polytric model of index one.
function lambdas = equi_polytrope_mass(N)

mr = @(ri) ((sin(pi*ri)-pi*ri.*cos(pi*ri))/pi);

lambdas = zeros(1, N);
for i = 1:N;
    lambdas(i) = fzero(@(ri) (mr(ri)-i/N), i/N);
end

lambdas = flip(lambdas);

end