function options = cmsset(varargin)
%CMSSET Create options structure used by CMSPlanet class methods.
%   OPTIONS = CMSSET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an options
%   structure OPTIONS in which the named properties have the specified values.
%   Any unspecified properties have default values. Case is ignored for property
%   names and unique partials are allowed.
%
%   CMSSET with no input or output arguments displays all property names and
%   their possible values.
%
%KNOWN PROPERTIES
%
%dJtol - Convergence tolerance for gravity coefficients [ positive real {1e-8} ]
%drhotol - Convergence tolerance for density adjustment [ positive real {1e-6} ]
%MaxIterBar - Number of iterations allowed for relaxation to barotrope [ positive integer {60} ]
%MaxIterHE - Number of iterations allowed for relaxation to equilibrium shape [ positive integer {60} ]
%xlayers - Solve shape functions on xlayers and spline the rest [ integer scalar or vector (-1 to disable) {-1} ]
%verbosity - Level of runtime messages [0 {1} 2 3 4]
%prerat - Precalculate powers of ratios of lambdas (trades memory for speed) [ {true} | false ]

% If no arguments print usage and return.
if (nargin == 0) && (nargout == 0)
    print_usage()
    return
end

% Define name-value pairs.
p = inputParser;
p.FunctionName = mfilename;

p.addParameter('dJtol',1e-8,@isposscalar)
p.addParameter('drhotol',1e-6,@isposscalar)
p.addParameter('MaxIterBar',60,@isposintscalar)
p.addParameter('MaxIterHE',60,@isposintscalar)
p.addParameter('xlayers',-1,@(x)validateattributes(x,{'numeric'},{'vector','integer'}))
p.addParameter('verbosity',1,@isnonnegintscalar)
p.addParameter('prerat',true,@islogicalscalar)

% undocumented or obsolete options
p.addParameter('prsmeth','linear'); % undocumented pressure interpolation method
p.addParameter('moimeth','midlayerz'); % undocumented moi integral method
p.addParameter('nangles',48,@isposintscalar)% #colatitudes used to define level surface
p.addParameter('kmax',30,@isposintscalar) % degree to carry out gravity mulitpole expansion
p.addParameter('TolX',1e-13,@isposscalar) % termination tolerance for root finding

% Parse name-value pairs and return.
p.parse(varargin{:})
options = p.Results;

end

function isposscalar(x)
validateattributes(x,{'numeric'},{'positive','scalar'})
end

function islogicalscalar(x)
validateattributes(x,{'logical'},{'scalar'})
end

function isposintscalar(x)
validateattributes(x,{'numeric'},{'positive','integer','scalar'})
end

function isnonnegintscalar(x)
validateattributes(x,{'numeric'},{'nonnegative','integer','scalar'})
end

function print_usage()
help(mfilename)
end
