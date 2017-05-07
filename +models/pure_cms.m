function cmp = pure_cms(lams, dels)
%PURE_CMS Directly set lambdas and deltas of underlying CMS object.
%   PURE_CMS(lams, dels) returns a CMSPlanet object whose cms field is initialized
%   with the properties lambdas and deltas equal to input arguments lams and dels,
%   respectively. The length of lams and dels (which of course must be vectors of
%   equal length) will determine the number of layers in the model.
%
%   Obviously this can be easily done in any script or in the Command Window but
%   some driver scripts require a handle to a parametrized model-building
%   function, so we supply one.

narginchk(2,2)

cmp = CMSPlanet(length(lams));
cmp.cms.lambdas = lams;
cmp.cms.deltas = dels;

end
