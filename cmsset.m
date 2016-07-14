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
%KNOWN PROPERTIES
%
%nlayers - Number of constant density layers [positive integer {512}]
%nangles - Number of colatitude points used to define level surfaces [positive integer {48}]
%kmax - Degree to carry out mulitpole expansion of gravity moments [positive even {12}]
%rcore - Core radius, normalized
%qrot - Dimensionless rotation parameter
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
p.addParameter('kmax',12,@isposintscalar)
p.addParameter('rcore',0.15,@isposnormalscalar)
p.addParameter('qrot',0,@isnonnegscalar)

% Parse name-value pairs and return.
p.parse(varargin{:})
options = p.Results;

end

function isposscalar(x)
validateattributes(x,{'numeric'},{'positive','scalar'})
end

function isnonnegscalar(x)
validateattributes(x,{'numeric'},{'nonnegative','scalar'})
end

function isposintscalar(x)
validateattributes(x,{'numeric'},{'positive','integer','scalar'})
end

function isposnormalscalar(x)
validateattributes(x,{'numeric'},{'nonnegative','scalar','<=',1})
end

function print_usage()
help(mfilename)
end
