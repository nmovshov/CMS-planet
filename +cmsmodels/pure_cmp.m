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
%
%   PURE_CMP(arho) is a single input variant using the n-by-2 array instead of
%   two vectors. This simplifies calls to some function-functions, like fminsearch
%   and mhsample. The first column is avec and the second column is rhovec.

narginchk(1,2)
if nargin == 1
    validateattributes(avec,{'numeric'},{'2d','ncols',2},1)
    rhovec = avec(:,2);
    avec = avec(:,1);
else
    validateattributes(avec,{'numeric'},{'vector'},1)
    validateattributes(rhovec,{'numeric'},{'vector'},2)
    assert(length(avec) == length(rhovec), 'length(avec) ~= length(rhovec)')
end

cmp = CMSPlanet(length(avec));
cmp.ai = avec;
drho = [rhovec(1); diff(rhovec)];
cmp.M = (4*pi/3)*(cmp.a0^3)*sum(drho.*cmp.cms.Vs);
cmp.rhoi = rhovec;

end
