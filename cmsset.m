function options = cmsset(varargin)
%CMSSET Create options structure for CMS-planet.
%   OPTIONS = CMSSET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an options
%   structure OPTIONS in which the named properties have the specified values.
%   Any unspecified properties have default values. Case is ignored for property
%   names.
%
%   CMSSET with no input arguments displays all property names and their
%   possible values.
%
%CMSSET PROPERTIES
%
%nlayers - Number of constant density layers [positive integer {512}]
%
%nangles - Number of colatitude points used to define level surfaces [positive integer {48}]
%
%nmoments - Degree to carry out mulitpole expansion of gravity moments [positive even {12}]
%
%   Note: defaults chosen to match Hubbard (2013) example.

% If no arguments print usage and return.
if (nargin == 0) && (nargout == 0)
    print_usage()
    return
end

% Define name-value pairs.
p = inputParser;
p.FunctionName = mfilename;

p.addParameter('nlayers',512,@isposintscalar)
p.addParameter('nangles',48,@isposintscalar)
p.addParameter('nmoments',12,@isposintscalar)

% Parse name-value pairs and return.
p.parse(varargin{:})
options = p.Results;

end

function isposintscalar(x)
validateattributes(x,{'numeric'},{'positive','integer','scalar'})
end

function print_usage()
fprintf('CMS-PLANET OPTIONS:\n');
fprintf('  nlayers: [positive scalar integer {512}]\n');
fprintf('  nangles: [positive scalar integer {48}]\n');
fprintf('  nmoments: [positive scalar even {12}]\n');
end
