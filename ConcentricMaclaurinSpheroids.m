classdef ConcentricMaclaurinSpheroids
    %CONCENTRICMACLAURINSPHEROIDS Implementation of CMS shape model.
    %   This class implements the iterative relaxation of concentric Maclaurin
    %   spheroids from starting (dimensionless) radii and densities to a
    %   self-consistent hydrostatic equilibrium shape, as explained in Hubbard
    %   (2013).
    
    %% Properties
    properties (GetAccess = public, SetAccess = private)
        opts    % holds CMS project-wide options
        lambdas % normalized layer equatorial radii
        deltas  % normalized density steps
        zetas   % normalized and scaled level-surface radii
        Js      % rescaled dimensionless gravity moments
    end
    
    %% The constructor
    methods
        function obj = ConcentricMaclaurinSpheroids(opts)
            if nargin == 0
                opts = cmsset();
            end
            obj.lambdas = linspace(1, opts.rcore, opts.nlayers)';
            obj.deltas = zeros(opts.nlayers, 1);
            obj.deltas(1) = 1;
            obj.zetas = ones(opts.nlayers, opts.nangles);
            obj.Js.tilde = zeros(opts.nlayers, opts.nmoments);
            obj.Js.tilde_prime = zeros(opts.nlayers, opts.nmoments);
            obj.Js.pprime = zeros(opts.nlayers, 1);
            obj.opts = opts;
        end
    end
    
    %% Ordinary methods
    methods
        function zeta_i = lvl_surf(obj,ilayer,mu)
            if ilayer == 0
                disp 'use eq. 50'
            else
                disp 'use eq. 51'
            end
        end
    end
    
    %% Static methods
    methods (Static)
        
    end % End of static methods
    
end % End of classdef

%% Class-related functions
function y = eq50(zeta0,mu,q,Jtwiddle,lambda,N,kmax)
% Equation 50 in Hubbard (2013)
x1 = 0;
for ii=1:N
    for kk=1:kmax
        k2 = 2*kk;
        x1 = x1 + Jtwiddle(ii,k2)*lambda(ii)^k2;
    end
end
x2 = 0;
for ii=1:N
    for kk=1:kmax
        k2 = 2*kk;
        x2 = x2 + Jtwiddle(ii,k2)*lambda(ii)^k2*zeta0^(-k2)*Pn(k2,mu);
    end
end
U0 = 1 + 0.5*q - x1;
U = 1/zeta0*(1 - x2) + 1/3*q*zeta0^2*(1 - Pn(2,mu));
y = U - U0;
end

function y = eq51(mu)
% Equation 51 in Hubbard (2013)

end

function y = Pn(n,x)
% Ordinary fully normalized Legendre polynomial
Pnm = legendre(n,x);
y = Pnm(1);
end
