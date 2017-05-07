function cmp = pure_cmp(avec, rhovec)
%PURE_CMP Directly set layer radii and densities of CMSPlanet.
%   PURE_CMP(avec, rhovec) returns a CMSPlanet object whose ai and rhoi properties
%   are initilized to input arguments avec and rhovec, respectively. The length of
%   avec and rhovec (which of course must be vectors of equal length) will
%   determine the number of layers in the model.
%
%   Obviously this can be easily done in any script or in the Command Window but
%   some driver scripts require a handle to a parametrized model-building
%   function, so we supply one.

narginchk(2,2)

cmp = CMSPlanet(length(avec));
cmp.ai = avec;
drho = [rhovec(1); diff(rhovec)];
cmp.M = (4*pi/3)*(cmp.a0^3)*sum(drho.*cmp.cms.Vs);
cmp.rhoi = rhovec;

end
