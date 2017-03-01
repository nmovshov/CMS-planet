%% A |CMSPlanet| quick start guide
% Example of interactive workflow with the |CMSPlanet| class. Note that most of
% the time it is better to use scripts to manage production runs. However the
% interactive workflow is useful for learning the basic capabilities of the class
% and is often more convenient for exploring initial ideas and results.

%% Create a |CMSPlanet| object
addpath(''); % path to your CMS-planet folder
N = 8; % Number of constant-density layers
cmp = CMSPlanet(N)

%%
% Notes:
%
% # The root CMS-planet folder must be on the MATLAB path
% # The "constructor" function returns a |CMSPlanet| object with many fields
% already initalized to sensible default values. Other fields are left empty and
% must be specified next.
% # The constructor can take many optional arguments. The one required argument is
% the number of layers to use. This number is now fixed.
% # For now the layers are equally spaced with normalized equatorial radii between
% |1| and |1/N|. This can be changed later.

%% Customize the |CMSPlanet| object
try
    si = zetUnits;
catch
    si = setFUnits;
end
cmp.M = 318*si.earth_mass;
cmp.a0 = 71490*si.km;
P_rot = 9*si.hour + 55*si.minute + 29.7*si.second;
w_rot = 2*pi/P_rot;
G = si.gravity;
cmp.qrot = w_rot^2*cmp.a0^3/(G*cmp.M);
cmp.name = 'dpl';
cmp.desc = 'a dummy planet';
cmp

%%
% Notes:
%
% # We set the total mass and outer equatorial radius of our planet. We also set
% the rotation parameter $q=\omega^2a_0^3/(GM)$.
% # You'll notice that many fields of |cmp| that were previously empty now have
% values. Many of the fields are self explanatory and you can also type |doc
% CMSPlanet| in the command window to see a list of fields and their definition.
% # Wondering about that business with |si=setUnits| and |si=setFUnits|? If you
% have the <http://github.com/nmovshov/physunits physunits toolbox> installed you
% can use it to benefit from automatic dimensional correctness and display. If you
% don't want to bother with it you can still use the function |setFUnits| to
% quickly get a struct variable with lots of predefined constants. But of course
% you don't have to use either.
% # Since we didn't specify a density distribution our planet is set by default
% with a constant density equal to the reference density $\rho_0 =
% \frac{3M}{4\pi{}a_0^3}$. All the layer densities, |rhoi|, are equal to each other
% and equal to |rho0|.
% # At this point our planet is not yet in hydrostatic equalibrium. It is still
% spherical, which is why the mean radius |s0| is equal to the equatorial radius
% |a0|, and why the mean density $\rho_s=\frac{3M}{4\pi{}s_0^3}$ is equal to the
% reference density $\rho_0$. You may have also noticed that the normalized moment
% of inertia, |NMoI|, is equal to |0.4| as expected for a homogenous sphere.

%% Relax to hydrostatic equilibrium
cmp.opts.verbosity = 2; % Control the verbosity of output, 0--4.
cmp.opts.dJtol = 1e-7; % For the demo we use a pretty lax tolerance.
cmp.relax_to_HE;
cmp.plot_equipotential_surfaces;
cmp

%%
% Notes:
%
% # The homogeneous planet is now in hydrostatic equilibrium and has aquired a
% noticeably flattened shape.
% # Many relevant fields have automatically updated their values. The mean radius
% is calculated by integrating the level surface of the outermost spheroid. The
% mean density, the gravity coefficients vector |Js|, the polar radius |b0|, the
% flattening |f0|, etc.
% # The pressure at the top of each layer, |Pi|, is calculated by integrating the
% hydrostatic equation. The pressure in the middle of each layer is in |P_mid|.
% The pressure in the center of the planet is in |P_c|.
% # Notice how the calculated mass |M_calc| is quite different from the mass we
% specified!

%% Match the specified total mass by adjusting the layer densities
% Because we fixed the layer densities and then changed their shapes to conform to
% equipotential surfaces the mass of the planet as calculated by a volume integral
% over all space is not the same as the value we specified initially. Usually we
% want the total mass to match a known value. We can do this by multiplying all
% densities by a constant factor. Doing this will adjust the total mass _without
% changing the gravity coefficients_. The convenience method |renormalize_density|
% will do this for us and remember the constant used in the public field
% |betanorm| in case we need it later.
cmp.renormalize_density;
cmp.betanorm

%% Specifying a density profile
% So far we worked with a constant density planet which is not that interesting.
% Let's now set a (slightly) more intersting density profile and see how it
% affects the equipotential surface and gravity moments. We can use any density
% profile we like. We just set the |rhoi| vector with |N| values. A physical
% density profile should be monotonically increasing towards the center.
cmp.rhoi = cmp.rho0*linspace(0, 3, cmp.nlayers);
cmp.opts.verbosity = 0;
cmp.relax_to_HE;
cmp.renormalize_density;
cmp.plot_barotrope;
cmp.plot_rho_of_r;
cmp

%% Converging to a desired barotrope
% We can ask |CMSPlanet| to find a density profile that will relax to a
% hydrosytatic configuration closely matching a desired pressure-density relation.
% Instead of specifying a density profile we specify a desired barotrope and save
% it to the field |eos|. Then we run the method |relax_to_barotrope()|. (Note that
% this process is more time consuming than simply relaxing to HE, because the
% density profile is constantly changing while the CMS algorithm is trying to find
% the equipotential surfaces!) We need an object derived from the abstract class
% |barotropes.Barotrope|.

help barotropes % It's easy to implement more barotropes
cmp.rhoi = cmp.rho0*ones(cmp.nlayers,1); % start with a constant density
cmp.eos = barotropes.Polytrope(2e5, 1);
cmp.opts.verbosity = 2;
cmp.opts.dBtol = 1e-7; % tolerance for density profile changes
cmp.relax_to_barotrope;
cmp.plot_barotrope('showinput', true, 'showscaledinput', true);
cmp.plot_rho_of_r;

%%
% Notes:
%
% # Before |relax_to_barotrope()| finished it called |renormalize_densities| so
% the mass of our planet matches the desired total mass already.
% # The pressure-density relation in our planet is not exactly the one we asked
% for but a closely related one. Since we had to scale the densities by a constant
% factor $\beta$ the resulting pressure-density relation changed. If we requested
% a barotrope of the form $P = f(\rho)$ we converge instead to the barotrope $P =
% \beta f(\beta^{-1}\rho).$ It's not always true that $\beta\approx{1}$, that
% depends on the match between the requested barotrope and the specified total
% mass.
