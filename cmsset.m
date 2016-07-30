function options = cmsset(varargin)
%CMSSET Create options structure used by CMS class methods.
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
%nlayers - Number of constant density layers [ positive integer {512} ]
%nangles - Number of colatitude points used to define level surfaces [ positive integer {48} ]
%kmax - Degree to carry out mulitpole expansion of gravity moments [ positive even {12} ]
%dJtol - Convergence tolerance for gravity moments [ positive real {1e-12} ]
%MaxIter - Maximum number of iterations allowed [ positive integer {40} ]
%rcore - Core radius, normalized [ {0.15} ]
%qrot - Dimensionless rotation parameter [ {0} ]
%J_integration_method - Choice of integration algorithm to compute J moments [ 'adaptive' | {'gauss'} ]
%zetas_in_J_integrals - How to obtain values of zeta inside J integrals [ {'rootfind'} | 'interp' ]
%verbosity - Level of runtime messages [0 {1} 2]
%IntTol - Relative tolerance for adaptive integrals [ positive real {1e-9} ]
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
p.addParameter('dJtol',1e-12,@isposscalar)
p.addParameter('MaxIter',40,@isposintscalar)
p.addParameter('rcore',0.15,@isposnormalscalar)
p.addParameter('qrot',0,@isnonnegscalar)
p.addParameter('verbosity',1,@isnonnegintscalar)
p.addParameter('J_integration_method','gauss',@isvalidintmethod)
p.addParameter('zetas_in_J_integrals','rootfind',@isvalidzetasmethod)
p.addParameter('IntTol',1e-9,@isposscalar)

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

function isnonnegintscalar(x)
validateattributes(x,{'numeric'},{'nonnegative','integer','scalar'})
end

function isposnormalscalar(x)
validateattributes(x,{'numeric'},{'nonnegative','scalar','<=',1})
end

function isvalidintmethod(x)
validateattributes(x,{'char'},{'row'})
if ~any(strcmpi(x,{'adaptive','gauss'}))
    error('Integration method %s not implemented.',x)
end
end

function isvalidzetasmethod(x)
validateattributes(x,{'char'},{'row'})
if ~any(strcmpi(x,{'rootfind','interp'}))
    error('Zetas in J integrals method %s not implemented.',x)
end
end

function print_usage()
help(mfilename)
end
