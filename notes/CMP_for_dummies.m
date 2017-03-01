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
