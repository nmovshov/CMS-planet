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
%
%   PURE_CMS(lamdel) is a single input variant using the n-by-2 array instead of
%   two vectors. This simplifies calls to some function-functions, like fminsearch
%   and mhsample. The first column is lambdas and the second column is deltas.

narginchk(1,2)
if nargin == 1
    validateattributes(lams,{'numeric'},{'2d','ncols',2},1)
    dels = lams(:,2);
    lams = lams(:,1);
else
    validateattributes(lams,{'numeric'},{'vector'},1)
    validateattributes(dels,{'numeric'},{'vector'},2)
    assert(length(lams) == length(dels), 'length(lams) ~= length(dels)')
end

cmp = CMSPlanet(length(lams));
cmp.cms.lambdas = lams;
cmp.cms.deltas = dels;

end
