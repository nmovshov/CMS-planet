function plot_cms(lambdas, deltas, zetas, mus)
% Visualize a CMS object

% Require R2016a to use the amazing polarplot features
if verLessThan('matlab','9')
    error('CMS plotting requires R2016a or later')
end

% Prepare colatitudes for polar plot
mu = [1, fliplr(mus), 0];
th = acos(mu);            % 0 to pi/2
th = [th, fliplr(pi-th)]; % 0 to pi
th = [th, (pi + th)];     % 0 to 2pi

% Prepare polar axes
figure;
pax = polaraxes;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
pax.ThetaAxisUnits = 'rad';
hold(pax, 'on')

% Plot level surfaces colored by layer density
cmap = parula;
rho = cumsum(deltas);
romin = min(rho); romax = max(rho);
lh = gobjects(size(lambdas));
bs = nan(size(lambdas));
for k=1:length(lh)
    xi = zetas(k,:)*lambdas(k);
    xi = [bs(k), fliplr(xi), lambdas(k)];
    xi = [xi, fliplr(xi)]; %#ok<AGROW>
    xi = [xi, fliplr(xi)]; %#ok<AGROW>
    lh(k) = polarplot(pax, th, xi);
    lh(k).Tag = 'equisurface';
    if (rho(k) <= romin)
        ci = 1;
    elseif (rho(k) >= romax)
        ci = length(cmap);
    else
        ci = fix((rho(k) - romin)/(romax - romin)*length(cmap)) + 1;
    end
    lh(k).Color = cmap(ci,:);
end

% Make outer surface more distinct
lh(1).LineWidth = 2;
lh(1).Color = 'k';

% Show grid lines above contours
pax.Layer = 'top';

end
