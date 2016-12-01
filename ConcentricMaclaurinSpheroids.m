classdef ConcentricMaclaurinSpheroids < handle
    %CONCENTRICMACLAURINSPHEROIDS Implementation of CMS shape model.
    %   This class implements the iterative relaxation of concentric Maclaurin
    %   spheroids from starting (dimensionless) radii and densities to a
    %   self-consistent hydrostatic equilibrium shape, as explained in Hubbard
    %   (2013).
    
    %% Properties
    properties (Access = public)
        opts    % holds CMS user configurable options
        lambdas % normalized layer equatorial radii
        deltas  % normalized density steps
    end
    properties (SetAccess = private)
        mus     % colatitude cosines
        zetas   % normalized and scaled level-surface radii
        Js      % rescaled, dimensionless, layer gravity moments
    end
    properties (Dependent)
        Vpu     % gravitational potential on fixed angles in planetary units
        Qpu     % rotation potential on fixed angles in planetary units
        Upu     % total potential on fixed angles in planetary units
        equiUpu % total potential at pole of (assumed) equipotential surfaces
        NMoI    % normalized moment of inertia
    end
    properties (Access = private)
        N            % real nlayers
        Pnmu         % values of Legendre polynomials at fixed colatitudes
        Pnzero       % values of Legendre polynomials at equator
        Pnone        % values of Legendre polynomials at pole
        gws          % weight factors for Gauss integration (correspond to mus)
        zeta1s       % normalized rescaled level-surface polar radii
        lamratpow    % stores values of lamdas-ratios-powers for quick retrieval
        realVpu      % stores values of Vpu for quick retrieval
        realequiU    % stores values of equiUpu for quick access
        os           % optimset struct for use by fzero
        cooked       % flag indicating obj.relax() was run
        fullyCooked  % flag indicating successful convergence
        realVpuMod   % flag triggering recalculation of realVpu
        realequiUMod % flag triggering recalculation of realequiU
    end
    properties (Dependent) % Convenience names
        nlayers % number of layers
        qrot    % dimensionless rotation parameter
        as      % normalized equatorial radii (another name for lambdas)
        bs      % normalized polar radii
        fs      % layer flattening, a.k.a, oblateness: (a - b)/a
        ars     % layer aspect ratio (b/a)
        ss      % layer mean radius normalized to a0
        Vs      % layer volume normalized to 4pi/3 a0^3
    end
    
    %% The constructor
    methods
        function obj = ConcentricMaclaurinSpheroids(nlay,varargin)
            %CONCENTRICMACLAURINSPHEROIDS Class constructor.
            
            % The constructor only dispatches to InitCMS().
            if nargin == 0
                error(['Required argument missing: specify number of layers',...
                    ' as first input argument.'])
            end
            validateattributes(nlay,...
                {'numeric'},{'positive','integer','scalar'},...
                '','nlayers',1)
            warning off CMS:obsolete
            op = cmsset(varargin{:});
            obj.opts = op; % (calls set.opts which calls cmsset(op) again)
            warning on CMS:obsolete
            
            obj.InitCMS(nlay, op);
        end
    end % End of constructor block
    
    %% Public methods
    methods (Access = public)
        function ET = relax(obj)
            % Iterate calculation of gravitational moments until converged.
            
            t_rlx = tic;
            
            % Optional communication
            verb = obj.opts.verbosity;
            if (verb > 0)
                fprintf('  ')
                fprintf('Relaxing CMS to self-consistent level surfaces...\n\n')
            end
            if (verb > 2)
                try
                    sbj = ['CMS.relax() started on ', getenv('computername')];
                    sendmail(obj.opts.email,sbj)
                catch
                end
            end
            
            % Main loop
            dJ = Inf;
            iter = 1;
            while (dJ > obj.opts.dJtol) && (iter <= obj.opts.MaxIterHE)
                t_pass = tic;
                if (verb > 0)
                    fprintf('    ')
                    fprintf('Pass %d (of max %d)...\n', iter,...
                        obj.opts.MaxIterHE)
                end
                obj.update_zetas;
                dJ = obj.update_Js();
                if (verb > 0)
                    fprintf('    ')
                    fprintf('Pass %d (of max %d)...done. (%g sec.)\n',...
                        iter, obj.opts.MaxIterHE, toc(t_pass))
                    if (verb > 1)
                        fprintf('    ')
                        fprintf('dJ = %g; required tolerance = %g.\n\n',...
                            dJ, obj.opts.dJtol)
                    else
                        fprintf('\n')
                    end
                end
                if (verb > 3)
                    try
                        sbj = ['CMS.relax() on ', getenv('computername')];
                        msg{1} = sprintf(...
                            'Pass %d (of max %d)...done. (%g sec)',...
                            iter, obj.opts.MaxIterHE, toc(t_pass));
                        msg{2} = sprintf(...
                            'dJ = %g; required tolerance = %g.',...
                            dJ, obj.opts.dJtol);
                        sendmail(obj.opts.email,sbj,msg)
                    catch
                    end
                end
                iter = iter + 1;
            end
            
            % Flags and maybe warnings
            obj.cooked = true;
            if (dJ < obj.opts.dJtol)
                obj.fullyCooked = true;
            else
                if (verb > 0)
                    msg = ['Object may not have fully converged.\n',...
                        'Try increasing the convergence tolerance ',...
                        '(%s.opts.dJtol) and/or iteration limit ',...
                        '(%s.opts.MaxIterHE).\n'];
                    warning off backtrace
                    warning(msg, inputname(1), inputname(1))
                    warning on backtrace
                end
            end
            
            % Calculate polar radii
            for ii=1:obj.nlayers
                obj.zeta1s(ii) = obj.zeta_j_of_mu(ii, 1);
            end
            
            ET = toc(t_rlx);
            
            % Optional communication
            if (verb > 0)
                msg = 'Relaxing CMS to self-consistent level surfaces...done.';
                fprintf('  ')
                fprintf([msg, '\n'])
                try
                    fprintf('  ')
                    fprintf('Total elapsed time %s\n',lower(seconds2human(ET)))
                catch
                    fprintf('  ')
                    fprintf('Total elapsed time %g sec.\n', ET)
                end
            end
            if (verb > 2)
                try
                    sbj = ['CMS.relax() finished on ', getenv('computername')];
                    sendmail(obj.opts.email,sbj)
                catch
                end
            end
        end
        
        function update_zetas(obj)
            % Update level surfaces using current value of Js.
            
            % Precompute powers of ratios of lambdas
            if isempty(obj.lamratpow)
                obj.lamratpow = nan(obj.opts.kmax+2, obj.N, obj.N);
                for ii=1:obj.N
                    for jj=1:obj.N
                        for kk=1:obj.opts.kmax+2
                            obj.lamratpow(kk,ii,jj) =...
                                (obj.lambdas(ii)/obj.lambdas(jj))^(kk-1);
                        end
                    end
                end
            end
            
            % Optional communication
            verb = obj.opts.verbosity;
            if (verb > 0)
                t_z_pass = tic;
                fprintf('    Updating zetas....')
            end
            
            % Loop over layers (outer) and colatitudes (inner)
            obj.os = optimset('TolX',obj.opts.TolX);
            nbLats = obj.opts.nangles;
            newzetas = NaN(size(obj.zetas));
            parfor ii=1:obj.nlayers
                for alfa=1:nbLats
                    newzetas(ii,alfa) = obj.zeta_j_of_alfa(ii,alfa); %#ok<PFBNS>
                end
            end
            obj.zetas = newzetas;
            
            % Optional communication
            if (verb > 1)
                t_z_pass = toc(t_z_pass);
                fprintf('done. (%g sec)\n', t_z_pass)
            elseif (verb > 0)
                fprintf('done.\n')
            end
            
            % Flag Vpu and equiUpu re-calc
            obj.realVpuMod = true;
            obj.realequiUMod = true;
        end
        
        function dJ = update_Js(obj)
            % Single-pass update of gravitational moments (dispatcher).
            
            % Optional communication
            verb = obj.opts.verbosity;
            if (verb > 0)
                t_J_pass = tic;
                fprintf('    Updating Js....')
            end
            
            % Dispatch based on chosen integration method.
            switch lower(obj.opts.J_integration_method)
                case 'gauss'
                    dJ = obj.update_Js_gauss();
                case 'adaptive'
                    dJ = obj.update_Js_adaptive();
                otherwise
                    error('Unknown integration method: %s.',...
                           obj.opts.J_integration_method)
            end
            obj.Js.tpprime = obj.Js.pprime.*obj.lambdas.^3;
            obj.Js.Jn = obj.Jn();
            
            % Optional communication
            if (verb > 1)
                t_J_pass = toc(t_J_pass);
                fprintf('done. (%g sec)\n', t_J_pass)
            elseif (verb > 0)
                fprintf('done.\n')
            end
            
            % Flag Vpu and equiUpu re-calc
            obj.realVpuMod = true;
            obj.realequiUMod = true;
        end
        
        function TF = validate(obj)
            % Run a few simple sanity checks, warn if fail.
            
            TF = true;
            warning('off','backtrace')
            
            % Warn but don't fail uncooked objects
            if (~obj.cooked)
                warning('Uncooked object, try running %s.relax()',inputname(1))
            elseif (~obj.fullyCooked)
                msg = ['Object may not be fully converged. ',...
                       'Try running %s.relax() and increasing ',...
                       'the convergence tolerance (%s.opts.dJtol) ',...
                       'and/or iteration limit (%s.opts.MaxIterHE).'];
                warning(msg, inputname(1), inputname(1), inputname(1))
            end
            
            % Check for properly normalized zetas
            for j=1:obj.nlayers
                er = abs(1 - obj.zeta_j_of_mu(j,0));
                if er > 1e-14
                    warning('Layer %d normalized radius deviates by %g',j,er)
                    TF = false;
                end
            end
            
            % Check for properly normalized Js
            er = obj.Jn(0) + 1;
            if abs(er) > 1e-15
                warning('J0 + 1 = %g (should be near zero)',er)
                TF = false;
            end
            
            % Check for weak concentricity
            ind = find(diff(obj.bs) > 0);
            for k=1:length(ind)
                msg = 'Concentricity may be broken: polar radius b_%d > b_%d';
                warning(msg,ind(k) + 1, ind(k))
                TF = false;
            end
            
            % Check for strong concentricity
            ind = find((obj.bs(1:end-1) - obj.lambdas(2:end)) < 0);
            for k=1:length(ind)
                msg = 'Strong concentricity may be broken: b_%d < a_%d';
                warning(msg,ind(k), ind(k) + 1)
                TF = false;
            end
            
            % Return.
            warning('on','backtrace')
        end
        
        function V = Vext(obj,xi,mu,nmax)
            % Convenience method: return potential at external point (xi,mu).
            
            validateattributes(xi,{'numeric'},{'real','positive','scalar'},...
                '','xi',1)
            validateattributes(mu,{'numeric'},{'real','scalar',...
                '>=',-1,'<=',1},'','mu',2)
            assert(xi >= obj.xi_i_of_mu(1,mu),...
                'The point (xi,mu) = (%g,%g) is not external.',xi,mu);
            if nargin < 4
                nmax = obj.opts.kmax;
            else
                validateattributes(nmax,{'numeric'},...
                    {'positive','integer','scalar','<=',obj.opts.kmax},...
                    '','nmax',3)
            end
            tk = 0:2:nmax;
            p2k = arrayfun(@Pn,tk,mu*ones(size(tk)));
            V = -(1/xi)*sum(obj.Jn(tk).*xi.^-tk.*p2k);
        end
        
        function J = Jn(obj,n)
            % Convenience method: return external, not re-scaled J multipoles.
            
            if exist('n','var')
                validateattributes(n,{'numeric'},...
                                     {'nonnegative','even','<=',obj.opts.kmax})
            else
                n = 0:2:obj.opts.kmax;
            end
            for k=1:length(n)
                J(k) = dot(obj.Js.tilde(:,n(k)+1),obj.lambdas.^n(k));%#ok<AGROW>
            end
        end
        
        function xi = xi_i_of_mu(obj,ilayer,mu)
            % Normalized (but not rescaled) layer shape function.
            
            validateattributes(ilayer,{'numeric'},...
                {'positive','integer','scalar','<=',length(obj.lambdas)},...
                '','ilayer',1)
            validateattributes(mu,{'numeric'},...
                {'real','scalar','>=',0,'<=',1},...
                '','mu',2)
            
            xi = obj.lambdas(ilayer)*obj.zeta_j_of_mu(ilayer, mu);
        end
        
        function ah = plot(obj)
            % Visualize a CMS object, return axes handle.
            
            % Require R2016a to use the amazing polarplot features
            if verLessThan('matlab','9')
                error('CMS plotting requires R2016a or later')
            end
                        
            % Prepare colatitudes for polar plot
            mu = [1, fliplr(obj.mus), 0];
            th = acos(mu);            % 0 to pi/2
            th = [th, fliplr(pi-th)]; % 0 to pi
            th = [th, (pi + th)];     % 0 to 2pi
                        
            % Prepare polar axes
            figure;
            pax = polaraxes;
            pax.ThetaZeroLocation = 'top';
            pax.ThetaDir = 'clockwise';
            pax.ThetaAxisUnits = 'rad';
            hold(pax, 'on')
            
            % Plot level surfaces colored by layer density
            cmap = parula;
            rho = cumsum(obj.deltas);
            romin = min(rho); romax = max(rho);
            lh = gobjects(size(obj.lambdas));
            for k=1:obj.nlayers
                xi = obj.zetas(k,:)*obj.lambdas(k);
                xi = [obj.bs(k), fliplr(xi), obj.lambdas(k)];
                xi = [xi, fliplr(xi)]; %#ok<AGROW>
                xi = [xi, fliplr(xi)]; %#ok<AGROW>
                lh(k) = polarplot(pax, th, xi);
                lh(k).Tag = 'equisurface';
                if (rho(k) <= romin)
                    ci = 1;
                elseif (rho(k) >= romax)
                    ci = length(cmap);
                else
                    ci = fix((rho(k) - romin)/(romax - romin)*length(cmap)) + 1;
                end
                lh(k).Color = cmap(ci,:);
            end
            
            % Make outer surface more distinct
            lh(1).LineWidth = 2;
            lh(1).Color = 'k';
            
            % Show grid lines above contours
            pax.Layer = 'top';
            
            % Fake a solid core color (xi and ci already set)
            for scl=0:0.05:1
                polarplot(pax, th, xi*scl, 'color', cmap(ci,:), 'tag', 'core');
            end
            
            % Return handle if requested
            if (nargout == 1), ah = pax; end            
        end
        
    end % End of public methods block
    
    %% Private methods
    methods (Access = private)
        function InitCMS(obj,nlay, op)
            % (Re)Initialize a CMS object.
            
            % The one required and immutable parameter
            obj.N = nlay;
            
            % Default layer setup is linear
            obj.lambdas = linspace(1, 1/nlay, nlay)';
            
            % Default deltas setup is constant density
            obj.deltas = zeros(nlay, 1);
            obj.deltas(1) = 1;
            
            % Default mu setup is almost never used!
            obj.mus = linspace(0,1,op.nangles); % will be modified by gauss
            
            % Default zetas setup is spherical
            obj.zetas = ones(nlay, op.nangles);
            obj.zeta1s = ones(nlay, 1);
            
            % Default Js setup is spherical
            obj.allocate_spherical_Js(nlay, op.kmax);
            
            % Get mus and weights for Gaussian quadrature
            if strcmpi(op.J_integration_method, 'gauss')
                [obj.mus, obj.gws] = gauleg(0, 1, op.nangles);
            end
            
            % Precompute Legendre polynomials for fixed colatitudes (gauss quad)
            for k = 0:op.kmax
                obj.Pnmu(k+1,1:op.nangles) = Pn(k, obj.mus);
                obj.Pnzero(k+1,1) = Pn(k, 0);
                obj.Pnone(k+1,1) = Pn(k, 1);
            end
            
            % Set flags and counters
            obj.cooked = false;
            obj.fullyCooked = false;
            obj.realVpuMod = true;
            obj.realequiUMod = true;
        end
        
        function dJ = update_Js_gauss(obj)
            % Single-pass update of gravitational moments by Gaussian quad.
                        
            % Do common denominator in eqs. (40-43)
            denom = 0;
            for j=1:obj.nlayers
                fun = obj.zetas(j,:).^3;
                I = obj.gws*fun'; % gauss quad formula
                denom = denom + obj.deltas(j)*obj.lambdas(j)^3*I;
            end
            
            % Do J tilde, eq. (40)
            new_tilde = zeros(size(obj.Js.tilde));
            for ii=1:obj.nlayers
                for kk=0:obj.opts.kmax
                    if rem(kk, 2), continue, end
                    fun = obj.Pnmu(kk+1,:).*obj.zetas(ii,:).^(kk+3);
                    I = obj.gws*fun'; % gauss quad formula
                    new_tilde(ii,kk+1) = -(3/(kk + 3))*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                end
            end
            
            % Do J tilde prime, eqs. (41,42)
            new_tprime = zeros(size(obj.Js.tilde_prime));
            for ii=1:obj.nlayers
                for kk=0:obj.opts.kmax
                    if rem(kk, 2), continue, end
                    if kk == 2
                        % eq. (42)
                        fun = obj.Pnmu(3,:).*log(obj.zetas(ii,:));
                        I = obj.gws*fun'; % gauss quad formula
                        new_tprime(ii,kk+1) = -3*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                    else
                        % eq. (41)
                        fun = obj.Pnmu(kk+1,:).*obj.zetas(ii,:).^(2 - kk);
                        I = obj.gws*fun'; % gauss quad formula
                        new_tprime(ii,kk+1) = -(3/(2 - kk))*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                    end
                end
            end
            
            % Do J double prime, eq. (27)
            new_pprime = zeros(size(obj.Js.pprime));
            denom = 0;
            for j=1:obj.nlayers
                fun = obj.zetas(j,:).^3;
                I = obj.lambdas(j)^3*obj.gws*fun';
                denom = denom + obj.deltas(j)*I;
            end
            denom = 2*denom;
            for ii=1:obj.nlayers
                new_pprime(ii) = obj.deltas(ii)/denom;
            end
            
            % Return max change in Js and update Js in obj.
            dJ = max(abs(new_tilde(:) - obj.Js.tilde(:)));
            dJ = max(dJ, max(abs(new_tprime(:) - obj.Js.tilde_prime(:))));
            dJ = max(dJ, max(abs(new_pprime(:) - obj.Js.pprime(:))));
            
            obj.Js.tilde = new_tilde;
            obj.Js.tilde_prime = new_tprime;
            obj.Js.pprime = new_pprime;
        end
        
        function dJ = update_Js_adaptive(obj)
            % Single-pass update of gravitational moments using built-in integral.
            
            % Do common denominator in eqs. (40-43)
            denom = 0;
            for j=1:obj.nlayers
                fun = @(mu)obj.zeta_j_of_mus(j, mu).^3;
                I = integral(fun, 0, 1, 'RelTol', obj.opts.IntTol);
                denom = denom + obj.deltas(j)*obj.lambdas(j)^3*I;
            end
            
            % Do J tilde, eq. (40)
            new_tilde = zeros(size(obj.Js.tilde));
            for ii=1:obj.nlayers
                for kk=0:obj.opts.kmax
                    if rem(kk, 2), continue, end
                    fun = @(mu)Pn(kk,mu).*obj.zeta_j_of_mus(ii, mu).^(kk+3);
                    I = integral(fun, 0, 1, 'RelTol', obj.opts.IntTol);
                    new_tilde(ii,kk+1) = -(3/(kk + 3))*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                end
            end
            
            % Do J tilde prime, eqs. (41,42)
            new_tprime = zeros(size(obj.Js.tilde_prime));
            for ii=1:obj.nlayers
                for kk=0:obj.opts.kmax
                    if rem(kk, 2), continue, end
                    if kk == 2
                        % eq. (42)
                        fun = @(mu)Pn(2,mu).*log(obj.zeta_j_of_mus(ii, mu));
                        I = integral(fun, 0, 1, 'RelTol', obj.opts.IntTol);
                        new_tprime(ii,kk+1) = -3*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                    else
                        % eq. (41)
                        fun = @(mu)Pn(kk,mu).*obj.zeta_j_of_mus(ii, mu).^(2 - kk);
                        I = integral(fun, 0, 1, 'RelTol', obj.opts.IntTol);
                        new_tprime(ii,kk+1) = -(3/(2 - kk))*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                    end
                end
            end
            
            % Do J double prime, eq. (27)
            new_pprime = zeros(size(obj.Js.pprime));
            denom = 0;
            for j=1:obj.nlayers
                fun = @(mu)obj.zeta_j_of_mus(j, mu).^3;
                I = obj.lambdas(j)*integral(fun, 0, 1, 'RelTol', obj.opts.IntTol);
                denom = denom + obj.deltas(j)*I;
            end
            denom = 2*denom;
            for ii=1:obj.nlayers
                new_pprime(ii) = obj.deltas(ii)/denom;
            end
            
            % Return max change in Js and update Js in obj.
            dJ = max(abs(new_tilde(:) - obj.Js.tilde(:)));
            dJ = max(dJ, max(abs(new_tprime(:) - obj.Js.tilde_prime(:))));
            dJ = max(dJ, max(abs(new_pprime(:) - obj.Js.pprime(:))));
            
            obj.Js.tilde = new_tilde;
            obj.Js.tilde_prime = new_tprime;
            obj.Js.pprime = new_pprime;
        end
                
        function y = zeta_j_of_mu(obj,jlayer,mu)
            % Find lvl surface of jth layer at colat mu.
            assert(jlayer > 0 && jlayer <= obj.nlayers)
            assert(mu >= 0 && mu <=1)
            if jlayer == 1
                fun = @(x)eq50(x,mu,obj.Js.tilde,obj.lambdas,obj.opts.qrot);
            else
                fun = @(x)eq51(x,jlayer,mu,obj.Js.tilde,obj.Js.tilde_prime,obj.Js.pprime,obj.lambdas,obj.opts.qrot);
            end
            if strcmpi(obj.opts.rootfinder,'fzero')
                %y = fzero(fun, [0.6, 1.02]);
                y = fzero(fun, 1, obj.os);
            else
                y = nzero(fun, 0.6, 1, obj.opts.TolX, eps);
            end
        end
        
        function y = zeta_j_of_alfa(obj,jlayer,alfa)
            % Find lvl surface of jth layer at colat mu(alfa).
            assert(jlayer > 0 && jlayer <= obj.nlayers)
            assert(alfa > 0 && alfa <= obj.opts.nangles)
            if jlayer == 1
                fun = @(x)obj.eq50_alfa(x,alfa);
            else
                fun = @(x)obj.eq51_alfa(x,jlayer,alfa);
            end
            if strcmpi(obj.opts.rootfinder,'fzero')
                %y = fzero(fun, [0.6, 1.02]);
                %y = fzero(fun, 1);
                y = fzero(fun, obj.zetas(jlayer, alfa), obj.os);
            else
                y = nzero(fun, 0.6, 1, obj.opts.TolX, eps);
            end
        end
        
        function y = eq50_alfa(obj,zeta0,alfa)
            % Equation 50 in Hubbard (2013) for fixed colatitudes.
            
            % Local variables
            nbLayers = obj.N;
            nbMoments = obj.opts.kmax;
            Jt = obj.Js.tilde;
            q = obj.opts.qrot;
            P0 = obj.Pnzero;
            Pmu = obj.Pnmu(:,alfa);
            rats = obj.lamratpow;
            
            % Double sum in eq. (47)
            x1 = 0;
            for ii=1:nbLayers
                for kk=2:2:nbMoments % (note ind shift, start ind, odd J=0)
                    x1 = x1 + Jt(ii,kk+1)*rats(kk+1,ii,1)*P0(kk+1);
                end
            end
            
            % Double sum in eq. (50)
            x2 = 0;
            for ii=1:nbLayers
                for kk=2:2:nbMoments % (note ind shift, start ind, odd J=0)
                    x2 = x2 + Jt(ii,kk+1)*rats(kk+1,ii,1)*zeta0^(-kk)*Pmu(kk+1);
                end
            end
            
            % And combine
            U0 = 1 + 0.5*q - x1;
            U = (1/zeta0)*(1 - x2) + 1/3*q*zeta0^2*(1 - Pmu(3));
            y = U - U0;
        end

        function y = eq51_alfa(obj,zeta_j,jj,alfa)
            % Equation 51 in Hubbard (2013) for fixed colatitudes.
            
            % Local variables
            nbLayers = obj.N;
            nbMoments = obj.opts.kmax;
            Jt = obj.Js.tilde;
            Jtp = obj.Js.tilde_prime;
            Jpp = obj.Js.pprime;
            lambda = obj.lambdas;
            q = obj.opts.qrot;
            P0 = obj.Pnzero;
            Pmu = obj.Pnmu(:,alfa);
            rats = obj.lamratpow;
            
            % Double sum, row 1
            x1 = 0;
            for ii=jj:nbLayers
                for kk=0:2:nbMoments % (note ind shift, start ind, odd J=0)
                    x1 = x1 + Jt(ii,kk+1)*rats(kk+1,ii,jj)*zeta_j^(-kk)*Pmu(kk+1);
                end
            end
            
            % Double sum, row 2
            x2 = 0;
            for ii=1:jj-1
                for kk=0:2:nbMoments
                    x2 = x2 + Jtp(ii,kk+1)*rats(kk+2,jj,ii)*zeta_j^(kk+1)*Pmu(kk+1);
                end
            end
            
            % Single sum, row 2
            x3 = 0;
            for ii=1:jj-1
                x3 = x3 + Jpp(ii)*rats(4,jj,1)*zeta_j^3;
            end
            
            % Double sum, row 3
            x4 = 0;
            for ii=jj:nbLayers
                for kk=0:2:nbMoments
                    x4 = x4 + Jt(ii,kk+1)*rats(kk+1,ii,jj)*P0(kk+1);
                end
            end
            
            % Double sum, row 4
            x5 = 0;
            for ii=1:jj-1
                for kk=0:2:nbMoments
                    x5 = x5 + Jtp(ii,kk+1)*rats(kk+2,jj,ii)*P0(kk+1);
                end
            end
            
            % Single sum, row 4
            x6 = 0;
            for ii=1:jj-1
                x6 = x6 + Jpp(ii)*rats(4,jj,1);
            end
            
            % And combine
            y = -(1/zeta_j)*(x1 + x2 + x3) + ...
                (1/3)*q*lambda(jj)^3*zeta_j^2*(1 - Pmu(3)) + ...
                (x4 + x5 + x6) - 0.5*q*lambda(jj)^3;
            
        end
        
        function y = zeta_j_of_mus(obj,jlayer,mus)
            % ML quad requires y=f(x) take and return vector...
            assert(jlayer > 0 && jlayer <= obj.nlayers)
            assert(all(mus >= 0) && all(mus <=1))
            y = nan(size(mus));
            for alfa=1:numel(mus)
                if jlayer == 1
                    fun = @(x)eq50(x,mus(alfa),obj.Js.tilde,obj.lambdas,obj.opts.qrot);
                else
                    fun = @(x)eq51(x,jlayer,mus(alfa),obj.Js.tilde,obj.Js.tilde_prime,obj.Js.pprime,obj.lambdas,obj.opts.qrot);
                end
                %y = fzero(fun, [0.6, 1.02]);
                y(alfa) = fzero(fun, 1);
            end
        end
        
        function allocate_spherical_Js(obj,nlay,nmom)
            obj.Js.tilde = zeros(nlay, (nmom+1));
            obj.Js.tilde_prime = zeros(nlay, (nmom+1));
            obj.Js.pprime = zeros(nlay, 1);
            obj.Js.tpprime = zeros(nlay, 1);
            obj.Js.Jn = zeros(1, nmom+1);
            den = sum(obj.deltas.*obj.lambdas.^3);
            obj.Js.tilde(:,1) = -(obj.deltas.*obj.lambdas.^3)/den;
            obj.Js.tilde_prime(:,1) = -1.5*(obj.deltas.*obj.lambdas.^3)/den;
            obj.Js.pprime(:) = 0.5*obj.deltas/sum(obj.deltas);
            obj.Js.tpprime = obj.Js.pprime.*(obj.lambdas.^3);
            obj.Js.Jn(1) = sum(obj.Js.tilde(:,1));
        end
    end % End of private methods block
    
    %% Access methods
    methods
        function set.opts(obj,val)
            
            % First filter through cmsset again.
            val = cmsset(val);
            
            % Don't change these horses mid-stream, it's just too messy.
            triggerFields = {'kmax','nangles'};
            if ~isempty(obj.opts) % for the call in the constructor
                for k=1:length(triggerFields)
                    if val.(triggerFields{k}) ~= obj.opts.(triggerFields{k})
                        msg = ['Changing %s of an existing object is just ',...
                            'too messy; create a new object instead.'];
                            error(msg,triggerFields{k})
                    end
                end
            end

            % Assign the new opts, unflag cooked, and return.
            obj.opts = val;
            obj.fullyCooked = false; %#ok<MCSUP>
        end
        
        function set.lambdas(obj,val)
            validateattributes(val,{'numeric'},{'vector','positive','finite',...
                '<=',1},'','lambdas')
            usval = sort((double(val)),'descend');
            if ~isequal(double(val),usval)
                warning off backtrace
                warning('lambdas array was sorted.')
                warning on backtrace
            end
            assert(numel(usval) == obj.nlayers,...
                'length(lambdas) = %g ~= nlayers = %g',...
                numel(usval),obj.nlayers) %#ok<MCSUP>
            assert(usval(1) == 1,['lambdas(1) = %g\n',...
                'radii must be normalized to outer layer'],usval(1))
            obj.lambdas = usval(:);
            obj.cooked = false; %#ok<MCSUP>
            obj.fullyCooked = false; %#ok<MCSUP>
            obj.realVpuMod = true; %#ok<MCSUP>
            obj.realequiUMod = true; %#ok<MCSUP>
        end
        
        function set.deltas(obj,val)
            validateattributes(val,{'numeric'},{'vector','real','finite'},...
                '','deltas')
            assert(numel(val) == obj.nlayers,...
                'length(deltas) = %g ~= nlayers = %g',...
                numel(val),obj.nlayers) %#ok<MCSUP>
            obj.deltas = double(val(:));
            obj.cooked = false; %#ok<MCSUP>
            obj.fullyCooked = false; %#ok<MCSUP>            
        end
        
        function val = get.nlayers(obj)
            val = obj.N;
        end
        
        function set.nlayers(~,~)
            msg = ['Changing number of layers in an existing CMS object ',...
                'makes no sense; create a new object instead.'];
            error(msg)
        end
        
        function val = get.qrot(obj)
            val = obj.opts.qrot;
        end
        
        function set.qrot(obj,val)
            obj.opts.qrot = val;
        end
        
        function val = get.bs(obj)
            val = obj.zeta1s.*obj.lambdas;
        end
        
        function val = get.as(obj)
            val = obj.lambdas;
        end
        
        function val = get.fs(obj)
            val = (obj.as - obj.bs)./obj.as;
        end
        
        function val = get.ars(obj)
            val = obj.bs./obj.as;
        end
        
        function val = get.Vs(obj)
            val = NaN(size(obj.lambdas));
            for j=1:obj.nlayers
                val(j) = obj.lambdas(j)^3*(obj.zetas(j,:).^3)*(obj.gws');
            end
        end
        
        function val = get.ss(obj)
            val = obj.Vs.^(1/3);
        end
        
        function val = get.Vpu(obj)
            if (~obj.realVpuMod)
                val = obj.realVpu;
            else
            val = NaN(size(obj.zetas));
            lam = obj.lambdas;
            zet = obj.zetas;
            til = obj.Js.tilde;
            tilp = obj.Js.tilde_prime;
            tilpp = obj.Js.tpprime;
            P2k = obj.Pnmu;
            n = (0:obj.opts.kmax)';
            for j=1:obj.nlayers
                for alfa=1:length(obj.mus)
                   val(j,alfa) = 0;
                   for i=j:obj.nlayers
                       val(j,alfa) = val(j,alfa) + ...
                           sum((lam(i)/lam(j)).^n.*til(i,:)'.*zet(j,alfa).^-n.*P2k(:,alfa));
                   end
                   for i=1:j-1
                       val(j,alfa) = val(j,alfa) + ...
                           sum((lam(j)/lam(i)).^(n+1).*tilp(i,:)'.*zet(j,alfa).^(n+1).*P2k(:,alfa));
                   end
                   for i=1:j-1
                       val(j,alfa) = val(j,alfa) + ...
                           (lam(j)/lam(i))^3*tilpp(i)*zet(j,alfa)^3;
                   end
                   val(j,alfa) = val(j,alfa)*(-1/(lam(j)*zet(j,alfa)));
                end
            end
            obj.realVpu = val;
            obj.realVpuMod = false;
            end
        end
        
        function val = get.Qpu(obj)
            val = NaN(size(obj.zetas));
            for j=1:obj.nlayers
                val(j,:) = (1/3)*obj.qrot*obj.lambdas(j)^2*obj.zetas(j,:).^2.*(1 - obj.Pnmu(3,:));
            end
        end
        
        function val = get.Upu(obj)
            val = obj.Vpu + obj.Qpu;
        end
        
        function val = get.equiUpu(obj)
            % Return potential on equipotential surface by sampling the pole.
            
            if (~obj.realequiUMod)
                val = obj.realequiU;
                return
            end
            
            val = zeros(obj.nlayers,1);
            lam = obj.lambdas;
            zet = obj.zeta1s;
            til = obj.Js.tilde;
            tilp = obj.Js.tilde_prime;
            tilpp = obj.Js.tpprime;
            P2k = obj.Pnone;
            n = (0:obj.opts.kmax)';
            for j=1:obj.nlayers
                V = 0;
                for i=j:obj.nlayers
                    V = V + sum((lam(i)/lam(j)).^n.*til(i,:)'.*zet(j).^-n.*P2k(:));
                end
                for i=1:j-1
                    V = V + sum((lam(j)/lam(i)).^(n+1).*tilp(i,:)'.*zet(j).^(n+1).*P2k(:));
                    V = V + (lam(j)/lam(i))^3*tilpp(i)*zet(j)^3;
                end
                V = V*(-1/(lam(j)*zet(j)));
                val(j) = V;
            end
            obj.realequiU = val;
            obj.realequiUMod = false;
        end
        
        function val = get.NMoI(obj)
            % C/Ma^2, see eq. 5 in Hubbard & Militzer 2016
            num = 0;
            den = 0;
            for j=1:obj.N
                num = num + obj.deltas(j)*...
                    ((obj.zetas(j,:)*obj.lambdas(j)).^5)*obj.gws';
                den = den + obj.deltas(j)*...
                    ((obj.zetas(j,:)*obj.lambdas(j)).^3)*obj.gws';
            end
            val = 2/5*num/den + 2/3*obj.Js.Jn(2);
        end
    end % End of access methods block
    
    %% Static methods
    methods (Static)
        function obj = loadobj(s)
            obj = ConcentricMaclaurinSpheroids(1);
            obj.N = s.N;
            obj.opts = s.opts;
            obj.lambdas = s.lambdas;
            obj.deltas = s.deltas;
            obj.mus = s.mus;
            obj.zetas = s.zetas;
            obj.zeta1s = s.zeta1s;
            obj.Js = s.Js;
            obj.Pnmu = s.Pnmu;
            obj.Pnzero = s.Pnzero;
            obj.Pnone = s.Pnone;
            obj.lamratpow = s.lamratpow;
            obj.gws = s.gws;
            obj.os = s.os;
            obj.cooked = s.cooked;
            obj.fullyCooked = s.fullyCooked;
            obj.realVpu = s.realVpu;
            obj.realequiU = s.realequiU;
            obj.realVpuMod = s.realVpuMod;
            obj.realequiUMod = s.realequiUMod;
        end
        
    end % End of static methods block
    
end % End of classdef

%% Class-related functions
function y = eq50(zeta0,mu,Jt,lambda,qrot)
% Equation 50 in Hubbard (2013)

% summation limits
nlayers = length(lambda);
kmax = size(Jt, 2) - 1;

% Legendre polynomials at 0 and mu
P0(kmax+1) = 0;
Pmu(kmax+1) = 0;
for k=0:kmax
    P0(k+1) = Pn(k, 0);
    Pmu(k+1) = Pn(k, mu);
end

% Double sum in eq. (47)
x1 = 0;
for ii=1:nlayers
    for kk=2:2:kmax % (note ind shift, start ind, odd J=0)
        x1 = x1 + Jt(ii,kk+1)*lambda(ii)^kk*P0(kk+1);
    end
end

% Double sum in eq. (50)
x2 = 0;
for ii=1:nlayers
    for kk=2:2:kmax % (note ind shift, start ind, odd J=0)
        x2 = x2 + Jt(ii,kk+1)*lambda(ii)^kk*zeta0^(-kk)*Pmu(kk+1);
    end
end

% And combine
U0 = 1 + 0.5*qrot - x1;
U = (1/zeta0)*(1 - x2) + 1/3*qrot*zeta0^2*(1 - Pmu(3));
y = U - U0;
end

function y = eq51(zeta_j,jj,mu,Jt,Jtp,Jpp,lambda,qrot)
% Equation 51 in Hubbard (2013)

% summation limits
nlayers = length(lambda);
kmax = size(Jt, 2) - 1;

% Legendre polynomials at 0 and mu
P0(kmax+1) = 0;
Pmu(kmax+1) = 0;
for k=0:kmax
    P0(k+1) = Pn(k, 0);
    Pmu(k+1) = Pn(k, mu);
end

% Double sum, row 1
x1 = 0;
for ii=jj:nlayers
    for kk=0:2:kmax % (note ind shift, start ind, odd J=0)
        x1 = x1 + Jt(ii,kk+1)*(lambda(ii)/lambda(jj))^kk*zeta_j^(-kk)*Pmu(kk+1);
    end
end

% Double sum, row 2
x2 = 0;
for ii=1:jj-1
    for kk=0:2:kmax
        x2 = x2 + Jtp(ii,kk+1)*(lambda(jj)/lambda(ii))^(kk+1)*zeta_j^(kk+1)*Pmu(kk+1);
    end
end

% Single sum, row 2
x3 = 0;
for ii=1:jj-1
    x3 = x3 + Jpp(ii)*lambda(jj)^3*zeta_j^3;
end

% Double sum, row 3
x4 = 0;
for ii=jj:nlayers
    for kk=0:2:kmax
        x4 = x4 + Jt(ii,kk+1)*(lambda(ii)/lambda(jj))^kk*P0(kk+1);
    end
end

% Double sum, row 4
x5 = 0;
for ii=1:jj-1
    for kk=0:2:kmax
        x5 = x5 + Jtp(ii,kk+1)*(lambda(jj)/lambda(ii))^(kk+1)*P0(kk+1);
    end
end

% Single sum, row 4
x6 = 0;
for ii=1:jj-1
    x6 = x6 + Jpp(ii)*lambda(jj)^3;
end

% And combine
y = -(1/zeta_j)*(x1 + x2 + x3) + ...
    (1/3)*qrot*lambda(jj)^3*zeta_j^2*(1 - Pmu(3)) + ...
    (x4 + x5 + x6) - 0.5*qrot*lambda(jj)^3;

end

function y = Pn(n,x)
% Fast implementation of ordinary Legendre polynomials of low degree.
switch n
    case 0
        y = ones(size(x));
    case 1
        y = x;
    case 2
        y = 0.5*(3*x.^2 - 1);
    case 3
        y = 0.5*(5*x.^3 - 3*x);
    case 4
        y = (1/8)*(35*x.^4 - 30*x.^2 + 3);
    case 5
        y = (1/8)*(63*x.^5 - 70*x.^3 + 15*x);
    case 6
        y = (1/16)*(231*x.^6 - 315*x.^4 + 105*x.^2 - 5);
    case 7
        y = (1/16)*(429*x.^7 - 693*x.^5 + 315*x.^3 - 35*x);
    case 8
        y = (1/128)*(6435*x.^8 - 12012*x.^6 + 6930*x.^4 - 1260*x.^2 + 35);
    case 9
        y = (1/128)*(12155*x.^9 - 25740*x.^7 + 18018*x.^5 - 4620*x.^3 + 315*x);
    case 10
        y = (1/256)*(46189*x.^10 - 109395*x.^8 + 90090*x.^6 - 30030*x.^4 + 3465*x.^2 - 63);
    case 11
        y = (1/256)*(88179*x.^11 - 230945*x.^9 + 218790*x.^7 - 90090*x.^5 + 15015*x.^3 - 693*x);
    case 12
        y = (1/1024)*(676039*x.^12 - 1939938*x.^10 + 2078505*x.^8 - 1021020*x.^6 + 225225*x.^4 - 18018*x.^2 + 231);
    otherwise
        assert(isvector(x))
        Pnm = legendre(n,x);
        y = Pnm(1,:);
        if ~isrow(x), y = y'; end
end
end

function [x,w] = gauleg(x1,x2,n)
%GAULEG Calculate abscissas and weights for Gauss-Legendre n-point quadrature.
%   [x,w] = GAULEG(x1,x2,n) returns the abscissas x and weights w that can be
%   used to evaluate the definite integral, I, of a function well approximated
%   by an (2n - 1) degree polynomial in the interval [x1,x2] using the
%   Gauss-Legendre formula:
%
%       I = sum(w.*f(x))
%
%   Algorithm
%     This function is based on the C++ implementation of a routine with the
%     same name in Numerical Recipes, 3rd Edition. But in several places I opt
%     for readability over performance, on the assumption that this function is
%     most likely to be called in a setup routine rather than in an inner-loop
%     computation.
%
%   Example
%     fun = @(x)sin(x);
%     [x,w] = gauleg(0,pi,6);
%     I_adaptive = integral(fun,0,pi)
%     I_gaussleg = sum(w.*fun(x))
%         
% Author: Naor Movshovitz (nmovshov at google dot com)
%         Earth and Planetary Sciences, UC Santa Cruz
%
% Reference: William H. Press, Saul A. Teukolsky, William T. Vetterling, and
% Brian P. Flannery. 2007. Numerical Recipes 3rd Edition: The Art of Scientific
% Computing (3 ed.). Cambridge University Press, New York, NY, USA.

% Input parsing and minimal assertions
narginchk(3,3)
nargoutchk(2,2)
validateattributes(x1,{'numeric'},{'scalar','finite','real'},1)
validateattributes(x2,{'numeric'},{'scalar','finite','real'},2)
validateattributes(n,{'numeric'},{'scalar','finite','integer','>=',2},3)
assert(x2 > x1, 'Interval must be positive.');

% Local variables
tol = 1e-14;
m = ceil(n/2);
xmid = (x1 + x2)/2;
dx = (x2 - x1);
x = NaN(1,n);
w = NaN(1,n);

% Main loop
for j=1:m
    % Get j-th root of Legendre polynomial Pn, along with Pn' value there.
    z = cos(pi*((j - 1) + 0.75)/(n + 0.5)); % initial guess for j-th root
    while true
        % Calculate Pn(z) and Pn-1(z) and Pn'(z)
        p = NaN(1,n+1);
        p(1) = 1;
        p(2) = z;
        for k=2:n
            pkm1 = p(k);
            pkm2 = p(k-1);
            pk = (1/k)*((2*k - 1)*z*pkm1 - (k - 1)*pkm2);
            p(k+1) = pk;
        end
        pn = p(end);
        pp = (n*p(end-1) - n*z*p(end))/(1 - z^2);
        
        % And now Newton's method (we are hopefully very near j-th root)
        oldz = z;
        z = z - pn/pp;
        if abs(z - oldz) < tol, break, end
    end
    
    % Now use j-th root to get 2 abscissas and weights
    x(j)     = xmid - z*dx/2; % Scaled abscissa left of center
    x(n+1-j) = xmid + z*dx/2; % Scaled abscissa right of center
    w(j)     = dx/((1 - z^2)*pp^2);
    w(n+1-j) = w(j);
end

% Verify and return
assert(all(isfinite(x)))
assert(all(isfinite(w)))
end

function x_root = nzero(fun, x1, x2, tolx, tolf)
% Bare-bones, optimized lion hunter in the interval [x1,x2].

f_l = fun(x1);
f_u = fun(x2);
if (abs(f_l) < tolf)
    x_root = x1;
    return
end
if (abs(f_u) < tolf)
    x_root = x2;
    return
end
if (f_l*f_u > 0), error('lion hunt error'), end
if (f_l < 0)
    x_l = x1;
    x_u = x2;
else
    x_l = x2;
    x_u = x1;
    f_u = f_l;
end
x_n = x_u;
f_n = f_u;
dx = (x_u - x_l)/2;
while (abs(dx) > tolx) && (abs(f_n) > tolf)
    dx = (x_u - x_l)/2;
    x_n = x_l + dx;
    f_n = fun(x_n);
    if (f_n < 0)
        x_l = x_n;
    else
        x_u = x_n;
    end
end
x_root = x_n;
end
