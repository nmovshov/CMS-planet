classdef ConcentricMaclaurinSpheroids < handle
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
        mus     % colatitude cosines
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
            obj.mus = linspace(0,1,opts.nangles);
            obj.zetas = ones(opts.nlayers, opts.nangles);
            obj.Js.tilde = zeros(opts.nlayers, (opts.kmax+1));
            obj.Js.tilde_prime = zeros(opts.nlayers, (opts.kmax+1));
            obj.Js.pprime = zeros(opts.nlayers, 1);
            obj.opts = opts;
        end
    end
    
    %% Ordinary methods
    methods
        function update_zetas(obj)
            % Update level surfaces using current value of Js.
            for ii=1:size(obj.zetas, 1)
                for alfa=1:length(obj.mus)
                    obj.zeta_j_of_mu(ii, obj.mus(alfa));
                end
            end
        end       
    end
    
    methods
        function y = zeta_j_of_mu(obj,ilayer,mu)
            % Find lvl surface of ith layer at colat mu.
            assert(ilayer > 0 && ilayer < obj.opts.nlayers)
            assert(mu >= 0 && mu <=1)
            if ilayer == 1
                fun = @(x)obj.eq50(x,mu,obj.Js.tilde,obj.lambdas);
                y = fzero(fun, [0.2, 1]);
            else
                disp 'use eq.51'
            end
        end
        function y = eq50(obj,zeta0,mu,Jt,lambda)
            % Equation 50 in Hubbard (2013)
            
            % Legendre polynomials at 0 and mu
            P0(obj.opts.kmax+1) = 0;
            Pmu(obj.opts.kmax+1) = 0;
            for k=0:obj.opts.kmax
                P0(k+1) = Pn(k, 0);
                Pmu(k+1) = Pn(k, mu);
            end
            
            % Double sum in eq. (47)
            x1 = 0;
            for ii=1:obj.opts.nlayers
                for kk=2:obj.opts.kmax % (note ind shift, start ind, odd J=0)
                    x1 = x1 + Jt(ii,kk+1)*lambda(ii)^kk*P0(kk+1);
                end
            end
            
            % Double sum in eq. (50)
            x2 = 0;
            for ii=1:obj.opts.nlayers
                for kk=2:obj.opts.kmax % (note ind shift, start ind, odd J=0)
                    x2 = x2 + Jt(ii,kk+1)*lambda(ii)^kk*zeta0^(-kk)*Pmu(kk+1);
                end
            end
            
            % And combine
            U0 = 1 + 0.5*obj.opts.qrot - x1;
            U = (1/zeta0)*(1 - x2) + 1/3*obj.opts.qrot*zeta0^2*(1 - Pmu(3));
            y = U - U0;
        end
    end
    
    %% Static methods
    methods (Static)
        
    end % End of static methods
    
end % End of classdef

%% Class-related functions

function y = eq51(mu)
% Equation 51 in Hubbard (2013)

end

function y = Pn(n,x)
% Ordinary fully normalized Legendre polynomial
Pnm = legendre(n,x);
y = Pnm(1);
end
