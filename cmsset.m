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
%nangles - Number of colatitude points used to define level surfaces [ positive integer {48} ]
%kmax - Degree to carry out mulitpole expansion of gravity moments [ positive even {30} ]
%dJtol - Convergence tolerance for gravity moments [ positive real {1e-10} ]
%dBtol - Convergence tolerance for barotrope adjustment [ positive real {1e-10} ]
%MaxIterHE - Maximum number of iterations allowed for relaxation to hydrostatic equilibrium [ positive integer {40} ]
%MaxIterBar - Maximum number of iterations allowed for relaxation to barotrope [ positive integer {40} ]
%J_integration_method - Choice of integration algorithm to compute J moments [ 'adaptive' | {'gauss'} ]
%zetas_in_J_integrals - How to obtain values of zeta inside J integrals [ {'rootfind'} | 'interp' ]
%verbosity - Level of runtime messages [0 {1} 2 3 4]
%IntTol - Relative tolerance for adaptive integrals [ positive real {1e-9} ]
%rootfinder - Algorithm choice for solving the zeta equations [ {'fzero'} | 'lionhunt' ]
%TolX - Termination tolerance for root finding algorithms [ positive real {1e-13} ]
%equipotential_squeeze - How to collapse U(i,mu) to U(i) on equipotential surfaces [ 'mean' | 'polar' | {'midlat'} ]

% If no arguments print usage and return.
if (nargin == 0) && (nargout == 0)
    print_usage()
    return
end

% Define name-value pairs.
p = inputParser;
p.FunctionName = mfilename;

p.addParameter('nangles',48,@isposintscalar)
p.addParameter('kmax',30,@isposintscalar)
p.addParameter('dJtol',1e-10,@isposscalar)
p.addParameter('dBtol',1e-10,@isposscalar)
p.addParameter('MaxIterHE',40,@isposintscalar)
p.addParameter('MaxIterBar',40,@isposintscalar)
p.addParameter('verbosity',1,@isnonnegintscalar)
p.addParameter('J_integration_method','gauss',@isvalidintmethod)
p.addParameter('zetas_in_J_integrals','rootfind',@isvalidzetasmethod)
p.addParameter('IntTol',1e-9,@isposscalar)
p.addParameter('rootfinder','fzero',@isvalidrootfinder)
p.addParameter('TolX',1e-13,@isposscalar)
p.addParameter('email','',@isvalidemail)
p.addParameter('equipotential_squeeze','midlat',@isvalidequiUsqueeze)

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
if strcmpi(x,'adaptive')
    warning off backtrace
    warning('CMS:obsolete','Using adaptive integration is not recommended.')
    warning on backtrace
end
end

function isvalidzetasmethod(x)
validateattributes(x,{'char'},{'row'})
if ~any(strcmpi(x,{'rootfind','interp'}))
    error('Zetas in J integrals method %s not implemented.',x)
end
end

function isvalidrootfinder(x)
validateattributes(x,{'char'},{'row'})
if ~any(strcmpi(x,{'fzero','lionhunt'}))
    error('Unknown root finding algorithm %s.',x)
end
end

function isvalidequiUsqueeze(x)
validateattributes(x,{'char'},{'row'})
if ~any(strcmpi(x,{'mean','polar','midlat'}))
    error('Unknown equipotential squeeze option %s.',x)
end
end

function TF = isvalidemail(x)
if isempty(x), TF = true; return, end
validateattributes(x,{'char'},{'row'})
validemail='[a-z_.1-9]+@[a-z_.1-9]+\.(com|net|edu)';
imatch=regexp(x,validemail);
if isempty(imatch) || ~isscalar(imatch) || imatch > 1
    error('Not a valid email address')
end
end

function print_usage()
help(mfilename)
end
