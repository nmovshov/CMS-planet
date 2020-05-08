classdef CMSPlanet < handle
    %CMSPLANET Interior model of rotating fluid planet.
    %   This class implements a model of a rotating fluid planet using Concentric
    %   Maclaurin Spheroids to calculate the hydrostatic equilibrium shape and
    %   resulting gravity field. A CMSPlanet object is defined by a given mass,
    %   equatorial radius, rotation parameter, and optionally a barotrope.
    
    %% Properties
    properties
        name   % model name
        mass   % reference mass
        radius % reference radius (equatorial!)
        P0     % reference pressure
        ai     % vector of equatorial radii (top down, a0=ai(1) is outer radius)
        rhoi   % vector of spheroid densities
        qrot   % rotation parameter, w^2a0^3/GM
        eos    % barotrope(s) (tip: help barotropes for options)
        bgeos  % optional background barotrope
        fgeos  % optional foreground barotrope
        opts   % holds user configurable options (tip: help cmsset)
    end
    properties (SetAccess = private)
        N      % convenience name for length(obj.si)
        CMS    % spheroid shape and other information returned by cms.m
        Js     % external gravity coefficients (returned by cms.m)
        betam  % mass renormalization factor (obj.mass/obj.M)
        alfar  % radius renormalization factor (obj.radius/obj.a0)
    end
    properties (Dependent)
        M      % calculated mass
        mi     % cumulative mass below ai
        a0     % calculated equatorial radius (another name for obj.ai(1))
        s0     % calculated mean radius
        si     % calculated spheroid mean radii
        rhobar % calculated mean density
        mrot   % rotation parameter referenced to s0
        J2     % convenience alias to obj.Js(2)
        J4     % convenience alias to obj.Js(3)
        J6     % convenience alias to obj.Js(4)
        J8     % convenience alias to obj.Js(5)
    end
    properties (GetAccess = private)
        G      % Gravitational constant
        u      % let's hold a units struct for convenience
    end
    
    %% A simple constructor
    methods
        function obj = CMSPlanet(varargin)
            % A simple constructor of CMSPlanet objects.
            % CMSPlanet('OPTION1', VALUE, 'OPTION2', VALUE2,...)
            
            % Populate options struct
            obj.opts = cmsset(varargin{:});
            
            % Init privates
            try
                obj.u = setFUnits;
                obj.G = obj.u.gravity;
            catch ME
                if isequal(ME.identifier,'MATLAB:UndefinedFunction')
                    error('Check your path, did you forget to run setws()?\n%s',ME.message)
                end
                rethrow(ME)
            end
        end
    end % End of constructor block
    
    %% Public methods
    methods (Access = public)
        function ET = relax_to_barotrope(obj)
            % Relax equilibrium shape, Js, and density simultaneously.
            
            % First some checks.
            if isempty(obj.eos)
                warning('CMSPLANET:assertion',...
                    'Set valid barotrope first (<obj>.eos = <barotrope>).')
                return
            end
            if numel(obj.renormalize()) < 2
                warning('CMSPLANET:assertion',...
                    'First set reference mass and equatorial radius.')
                return
            end
            if isempty(obj.qrot)
                warning('CMSPLANET:assertion',...
                    'First set rotation parameter (<obj>.qrot).')
                return
            end
            if isempty(obj.P0)
                warning('CMSPLANET:P0',...
                    'Setting reference pressure to one bar (<obj>.P0=1e5).')
                obj.P0 = 1*obj.u.bar;
            end
            
            % Optional communication
            verb = obj.opts.verbosity;
            if (verb > 0)
                fprintf('Relaxing to desired barotrope...\n\n')
            end
            
            % Ready, set,...
            mihe = obj.opts.MaxIterHE;
            obj.opts.MaxIterHE = 2;
            warning('off','CMS:maxiter')
            t_rlx = tic;
            
            % Main loop
            iter = 1;
            while (iter <= obj.opts.MaxIterBar)
                t_pass = tic;
                if (verb > 0)
                    fprintf('Baropass %d (of max %d)...\n',...
                        iter, obj.opts.MaxIterBar)
                end
                
                old_Js = obj.Js;
                old_Js(abs(old_Js) < 1e-12) = 0; % a hack for nonrotating planets
                if isempty(old_Js), old_Js = [-1,zeros(1,15)]; end
                old_ro = obj.rhoi;
                obj.relax_to_HE();
                obj.update_densities;
                dJs = abs((obj.Js - old_Js)./(old_Js));
                dJs = max(dJs(isfinite(dJs(1:6))));
                dro = obj.rhoi./old_ro;
                dro = var(dro(isfinite(dro)));
                
                if (verb > 0)
                    fprintf('Baropass %d (of max %d)...done. (%g sec)\n',...
                        iter, obj.opts.MaxIterBar, toc(t_pass))
                    fprintf('var(drho) = %g (%g required); dJ = %g (%g required).\n\n',...
                        dro, obj.opts.drhotol, dJs, obj.opts.dJtol)
                end
                
                % The stopping criterion is to satisfy both J and rho tolerance
                if (dro < obj.opts.drhotol) && dJs < obj.opts.dJtol
                    break
                end
                
                % end the main loop
                iter = iter + 1;
            end
            ET = toc(t_rlx);
            
            % Renormalize and record renorm factors
            renorms = obj.renormalize();
            obj.alfar = renorms(1);
            obj.betam = renorms(2);
            
            % Some clean up
            obj.opts.MaxIterHE = mihe;
            warning('on','CMS:maxiter')
            
            % Optional communication
            if (verb > 0)
                fprintf('Relaxing to desired barotrope...done.\n')
                fprintf('Total elapsed time %s\n',lower(seconds2human(ET)))
            end
        end
        
        function [ET, dJ] = relax_to_HE(obj)
            % Call cms to obtain equilibrium shape and gravity.
            
            if (obj.opts.verbosity > 1)
                fprintf('  Relaxing to hydrostatic equilibrium...\n')
            end
            
            t_rlx = tic;
            zvec = obj.ai/obj.ai(1);
            dvec = obj.rhoi/obj.rhoi(end);
            if isempty(obj.CMS), obj.CMS.JLike = struct(); end
            [obj.Js, obj.CMS] = cms(zvec, dvec, obj.qrot,...
                'tol', obj.opts.dJtol, 'maxiter', obj.opts.MaxIterHE,...
                'xlayers', obj.opts.xlayers, 'J0s', obj.CMS.JLike,...
                'prerat', obj.opts.prerat);
            ET = toc(t_rlx);
            dJ = obj.CMS.dJs;
            
            if (obj.opts.verbosity > 1)
                fprintf('  Relaxing to hydrostatic equilibrium...done.\n')
                fprintf('  Elapsed time %s\n',lower(seconds2human(ET)))
            end
        end
        
        function dro = update_densities(obj)
            % Set layer densities to match prescribed barotrope.
            
            if isempty(obj.eos)
                warning('Make sure input barotrope (<obj>.eos) is set.')
                return
            end
            
            t_rho = tic;
            verb = obj.opts.verbosity;
            if (verb > 1)
                fprintf('  Updating layer densities...')
            end
            P = obj.P_mid();
            if isscalar(obj.eos)
                newro = obj.eos.density(P);
            else
                newro = repmat(obj.rhoi(1), obj.N, 1);
                for k=1:length(newro)
                    newro(k) = obj.eos(k).density(P(k));
                end
            end
            dro = ((newro - obj.rhoi)./obj.rhoi);
            if (verb > 2)
                fprintf('done. (%g sec)\n', toc(t_rho))
            elseif (verb > 1)
                fprintf('done.\n')
            end
            obj.rhoi = newro;
        end
        
        function P_m = P_mid(obj)
            % Mid-layer pressure (by interpolation)
            try
                P = obj.Pi;
                r = obj.ai;
                P_m = NaN(size(P));
                P_m(1:end-1) = (P(1:end-1) + P(2:end))/2;
                P_m(end) = P(end) + ...
                    (P(end) - P(end-1))/(r(end-1) - r(end))*(r(end)/2);
            catch ME
                warning(ME.identifier, ...
                    'Pressure integration failed with message %s',ME.message)
                P_m = [];
            end
        end
        
        function P_c = P_center(obj)
            % Central pressure (by extrapolation)
            try
                P_c = spline(obj.CMS.lambdas, obj.Pi, 0);
            catch
                P_c = [];
            end
        end
        
        function obj = set_observables(obj, obs)
            % Copy physical properties from an +observables struct.
            obj.mass = obs.M;
            obj.radius = obs.a0;
            obj.qrot = obs.q;
            obj.P0 = obs.P0;
            try
                obj.bgeos = obs.bgeos;
                obj.fgeos = obs.fgeos;
            catch
            end
        end

        function obj = set_J_guess(obj, jlike)
            % Use with caution
            obj.CMS.JLike = jlike;
        end
        
        function ab = renormalize(obj)
            % Match input and calculated mass and equatorial radius.

            try
                a = obj.radius/obj.a0;
                obj.ai = obj.ai*a;
            catch
                a = [];
            end
            try
                b = obj.mass/obj.M;
                obj.rhoi = obj.rhoi*b;
            catch
                b = [];
            end
            ab = [a, b];
        end
        
        function obj = fix_radius(obj)
            % Resize planet to match equatorial radius to observed value.
            
            if isempty(obj.radius) || isempty(obj.a0) || isempty(obj.ai)
                warning('Missing information; no action.')
                return
            end
            obj.ai = obj.ai*obj.radius/obj.a0;
        end
        
        function obj = fix_mass(obj)
            % Rescale density to match planet mass to observed value.
            
            if isempty(obj.mass) || isempty(obj.rhoi)
                warning('Missing information; no action.')
                return
            end
            obj.rhoi = obj.rhoi*obj.mass/obj.M;
        end
    end % End of public methods block
    
    %% Visualizers
    methods (Access = public)
        function [ah, lh] = plot_rho_of_r(obj, varargin)
            % Plot rho(r) where r is equatorial radius.
            
            % Don't bother if uninitialized
            if isempty(obj.ai) || isempty(obj.rhoi)
                warning('Uninitialized object.')
                return
            end
            
            % Input parsing
            p = inputParser;
            p.addParameter('axes', [], @(x)isscalar(x) && isgraphics(x, 'axes'))
            p.addParameter('plottype', 'line', @(x)isrow(x) && ischar(x))
            p.parse(varargin{:})
            pr = p.Results;
            
            % Prepare the canvas
            if isempty(pr.axes)
                fh = figure;
                set(fh, 'defaultTextInterpreter', 'latex')
                set(fh, 'defaultLegendInterpreter', 'latex')
                ah = axes;
                hold(ah, 'on')
            else
                ah = pr.axes;
                hold(ah, 'on')
            end
            
            % Prepare the data
            x = [obj.ai/obj.a0; 0];
            y = [obj.rhoi; obj.rhoi(end)];
            
            % Plot the lines (density in 1000 kg/m^3)
            if isequal(lower(pr.plottype), 'stairs')
                lh = stairs(x, y/1000);
            elseif isequal(lower(pr.plottype), 'line')
                lh = line(x, y/1000);
            else
                error('Unknown value of parameter plottype.')
            end
            lh.LineWidth = 2;
            if isempty(pr.axes)
                lh.Color = [0.31, 0.31, 0.31];
            end
            if isempty(obj.name)
                lh.DisplayName = 'CMS model';
            else
                lh.DisplayName = obj.name;
            end
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.Box = 'on';
                xlabel('Level surface equatorial radius, $a/a_0$', 'fontsize', 12)
                ylabel('$\rho$ [1000 kg/m$^3$]', 'fontsize', 12)
            else
                xlim('auto')
            end
        end
        
        function [ah, lh] = plot_barotrope(obj, varargin)
            % Plot P(rho) of current model and optionally of input barotrope.
            
            % Don't bother if there is no pressure
            if isempty(obj.Pi)
                warning('Uninitialized object. Remember to set obj.P0?')
                return
            end
            
            % Input parsing
            p = inputParser;
            p.addParameter('axes', [],...
                @(x)isscalar(x) && isgraphics(x,'axes') && isvalid(x));
            p.addParameter('showinput', false,...
                @(x)isscalar(x) && islogical(x));
            p.addParameter('showscaledinput', false,...
                @(x)isscalar(x) && islogical(x));
            p.addParameter('includecore', false,...
                @(x)isscalar(x) && islogical(x));
            p.parse(varargin{:})
            pr = p.Results;
            
            % Prepare the canvas
            if isempty(pr.axes)
                fh = figure;
                set(fh, 'defaultTextInterpreter', 'latex')
                set(fh, 'defaultLegendInterpreter', 'latex')
                ah = axes;
                hold(ah, 'on')
            else
                ah = pr.axes;
                hold(ah, 'on')
            end
            
            % Prepare the data: model
            x_cms = obj.rhoi;
            y_cms = obj.P_mid();
            
            % Prepare the data: input
            if pr.showinput && ~isempty(obj.eos) && (range(x_cms) > 0)
                x_bar = linspace(min(x_cms), max(x_cms));
                if isscalar(obj.eos)
                    y_bar = double(obj.eos.pressure(x_bar));
                else
                    v = 1:length(unique(x_cms));
                    ind = interp1(unique(x_cms), v, x_bar, 'nearest', 'extrap');
                    y_bar = nan(size(x_bar));
                    for k=1:length(x_bar)
                        y_bar(k) = double(obj.eos(ind(k)).pressure(x_bar(k)));
                    end
                end
            else
                y_bar = NaN;
            end
            
            % Prepare the data: scaled input
            if pr.showscaledinput && ~isempty(obj.eos) && (range(x_cms) > 0)
                x_bar = linspace(min(x_cms), max(x_cms));
                bnorm = obj.betam; % the mass renorm factor
                anorm = obj.alfar; % the radius renorm factor
                if isempty(bnorm), bnorm = nan; end
                if isempty(anorm), anorm = nan; end
                if isscalar(obj.eos)
                    y_bar_scl = double(bnorm/anorm*obj.eos.pressure(x_bar/bnorm));
                else
                    v = 1:length(unique(x_cms));
                    ind = interp1(unique(x_cms), v, x_bar, 'nearest', 'extrap');
                    y_bar_scl = nan(size(x_bar));
                    for k=1:length(x_bar)
                        y_bar_scl(k) = double(...
                            bnorm/anorm*obj.eos(ind(k)).pressure(x_bar(k)/bnorm));
                    end
                end
            else
                y_bar_scl = NaN;
            end
            
            % Plot the lines (pressure in GPa)
            lh(1) = stairs(x_cms, y_cms/1e9);
            lh(1).LineWidth = 2;
            if isempty(pr.axes)
                lh(1).Color = [0.31, 0.31, 0.31];
            end
            if isempty(obj.name)
                lh(1).DisplayName = 'CMS model';
            else
                lh(1).DisplayName = obj.name;
            end
            
            if pr.showinput && any(isfinite(y_bar))
                lh(end+1) = line(x_bar, y_bar/1e9);
                lh(end).Color = 'r';
                lh(end).LineStyle = '--';
                lh(end).DisplayName = 'input barotrope';
            end
            
            if pr.showscaledinput && any(isfinite(y_bar_scl))
                lh(end+1) = line(x_bar, y_bar_scl/1e9);
                lh(end).Color = [0, 0.5, 0];
                lh(end).LineStyle = '--';
                lh(end).DisplayName = 'input barotrope ($\frac{\beta}{\alpha}$-scaled)';
            end
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.Box = 'on';
                if (range(x_cms) > 0)
                    xlim([min(x_cms),max(x_cms)])
                end
                xlabel('$\rho$ [kg/m$^3$]')
                ylabel('$P$ [GPa]')
            else
                xlim('auto')
            end
            
            % Legend
            legend(ah, 'off')
            gh = legend(ah, 'show','location','nw');
            gh.FontSize = 11;
        end
        
        function [ah, lh] = plot_spheroid_J_contributions(obj, n, varargin)
            % Plot relative weights of *spheroids* contribution to Js.
            
            % Input parsing
            if nargin < 2
                fprintf('Usage:\n\tCMP.plot_spheroid_J_contributions([2,4,...]).\n')
                return
            end
            validateattributes(n,{'numeric'},{'row','positive','integer','even'})
            p = inputParser;
            p.addParameter('axes',[],@(x)isscalar(x)&&isgraphics(x, 'axes'))
            p.addParameter('cumulative',false,@(x)isscalar(x)&&islogical(x))
            p.parse(varargin{:})
            pr = p.Results;
            
            % Prepare the data
            y = nan(obj.N, length(n));
            for k=1:length(n)
                cJi = cumsum(obj.CMS.JLike.fulltilde(:,n(k)+1).*obj.CMS.lambdas.^n(k));
                dJi = sdderiv(obj.CMS.lambdas, cJi);
                if pr.cumulative
                    y(:,k) = cJi/cJi(end);
                else
                    y(:,k) = abs(dJi)/max(abs(dJi));
                end
            end
            
            % Prepare the canvas
            if isempty(pr.axes)
                fh = figure;
                set(fh, 'defaultTextInterpreter', 'latex')
                set(fh, 'defaultLegendInterpreter', 'latex')
                ah = axes;
            else
                ah = pr.axes;
                axes(ah)
            end
            hold(ah, 'on')
            
            % Plot the lines
            lh = gobjects(1,length(n));
            for k=1:length(n)
                lh(k) = plot(obj.CMS.lambdas, y(:,k), 'LineWidth',2);
                lh(k).DisplayName = sprintf('$J_{%i}$',n(k));
            end
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.Box = 'on';
                xlabel('Spheroid normalized equatorial radius, $\lambda$', 'fontsize', 12)
                if pr.cumulative
                    ylabel('$|J_n(\lambda>\lambda_i)|$ [normalized]', 'fontsize', 12)
                else
                    ylabel('$|J_n(\lambda_i)|$ [normalized]', 'fontsize', 12)
                end
            end
            
            % Legend
            legend(ah, 'off')
            gh = legend(ah, flipud(lh));
            if pr.cumulative
                gh.Location = 'southwest';
            else
                gh.Location = 'northwest';
            end
            gh.FontSize = 11;
        end
        
        function [ah, lh] = plot_moi_contribution(obj, varargin)
            % Plot relative contribution to MoI by depth.
            
            p = inputParser;
            p.addParameter('axes',[],@(x)isscalar(x)&&isgraphics(x, 'axes'))
            p.addParameter('cumulative',false,@(x)isscalar(x)&&islogical(x))
            p.parse(varargin{:})
            pr = p.Results;
            
            % Prepare the data
            if pr.cumulative
                y = obj.NMoI('csum');
                y = y/y(end);
            else
                y = obj.NMoI('none');
                y = y/max(abs(y));
            end
            
            % Prepare the canvas
            if isempty(pr.axes)
                fh = figure;
                set(fh, 'defaultTextInterpreter', 'latex')
                set(fh, 'defaultLegendInterpreter', 'latex')
                ah = axes;
            else
                ah = pr.axes;
                axes(ah)
            end
            hold(ah, 'on')
            
            % Plot the lines
            lh = plot(obj.si/obj.s0, y, 'LineWidth',2);
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.Box = 'on';
                xlabel('Spheroid normalized equatorial radius, $\lambda$', 'fontsize', 12)
                if pr.cumulative
                    ylabel('$I(\lambda>\lambda_i)$ [normalized]', 'fontsize', 12)
                else
                    ylabel('$I(\lambda_i)$ [normalized]', 'fontsize', 12)
                end
            end
        end
    end % End of visulaizers block
    
    %% Reporters/exporters
    methods (Access = public)
        function T = report_card(obj, obs)
            % REPORT_CARD Table summary of model's vital statistics.
            
            % Minimal checks
            narginchk(1,2);
            try
                obj.J2;
            catch
                warning('Uncooked object.')
                return
            end
            
            % Basic table
            vitals = {'Mass [kg]'; 'J2'; 'J4'; 'J6'; 'J8'; 'NMoI'};
            CMP1 = [obj.M; obj.J2; obj.J4; obj.J6; obj.J8; obj.NMoI];
            CMP1 = double(CMP1);
            T = table(CMP1, 'RowNames', vitals);
            if ~isempty(obj.name)
                vname = matlab.lang.makeValidName(obj.name);
                T.Properties.VariableNames{'CMP1'} = vname;
            end
            if nargin == 1, return, end
            
            % Optionally compare with something
            try
                oM = obs.M;
            catch
                oM = NaN;
            end
            try
                oJ2 = obs.J2;
                oJ4 = obs.J4;
                oJ6 = obs.J6;
                oJ8 = obs.J8;
            catch
                oJ2 = NaN;
                oJ4 = NaN;
                oJ6 = NaN;
                oJ8 = NaN;
            end
            try
                oNMoI = obs.NMoI;
            catch
                oNMoI = NaN;
            end
            try
                oname = obs.name;
            catch
                oname = [];
            end
            OBS1 = [oM; oJ2; oJ4; oJ6; oJ8; oNMoI];
            OBS1 = double(OBS1);
            T = [T table(OBS1)];
            if ~isempty(oname)
                vname = matlab.lang.makeValidName(obs.name);
                try
                    T.Properties.VariableNames{'OBS1'} = vname;
                catch
                    T.Properties.VariableNames{'OBS1'} = ['x_',vname];
                end
            end
            DIFF = (CMP1 - OBS1)./CMP1;
            T = [T table(DIFF, 'VariableNames', {'frac_diff'})];
        end
        
        function s = to_struct(obj, rdc, keepjlike)
            % Convert object to static struct keeping only essential fields.
            
            if nargin < 2, rdc = 1; end % 0=none, 1=to double, 2=to single, 3=to scalars
            if nargin < 3, keepjlike = false; end % keep JLike e.g. to help re-relaxing
            
            s.name   = obj.name;
            s.M      = obj.M;
            s.s0     = obj.s0;
            s.a0     = obj.a0;
            s.M_Z    = obj.M_Z;
            s.rhobar = obj.rhobar;
            s.mrot   = obj.mrot;
            s.qrot   = obj.qrot;
            s.J2     = obj.J2;
            s.J4     = obj.J4;
            s.J6     = obj.J6;
            s.J8     = obj.J8;
            s.NMoI   = obj.NMoI;
            s.si     = obj.si;
            s.ai     = obj.ai;
            s.rhoi   = obj.rhoi;
            s.Pi     = obj.Pi;
            s.mi     = obj.mi;
            s.zi     = obj.zi;
            
            if rdc > 0
                s = structfun(@double, s, 'UniformOutput', false);
                s.name = obj.name;
            end
            if rdc > 1
                s = structfun(@single, s, 'UniformOutput', false);
                s.name = obj.name;
            end
            if rdc > 2
                s.si     = [];
                s.ai     = [];
                s.rhoi   = [];
                s.Pi     = [];
                s.mi     = [];
                s.zi     = [];
            end
            
            try
                s.eos = obj.eos.name;
            catch
                s.eos = '';
            end
            try
                s.bgeos = obj.bgeos.name;
            catch
                s.bgeos = '';
            end
            try
                s.fgeos = obj.fgeos.name;
            catch
                s.fgeos = '';
            end
            
            if keepjlike
                s.CMS.JLike = obj.CMS.JLike;
            else
                s.CMS.JLike = [];
            end
        end
        
        function T = to_table(obj)
            % Return a table of critical quantities.
            
            T = table;
            T.ai = double(obj.ai);
            T.si = double(obj.si);
            T.rhoi = double(obj.rhoi);
            T.Pi = double(obj.Pi);
            T.mi = double(obj.mi);
            if ~isempty(obj.zi)
                T.zi = double(obj.zi);
            end
        end
        
        function to_ascii(obj, fname)
            % Export the state of the model as ascii file.
            
            % File name
            if nargin < 2
                fprintf('Usage:\n\tcmp.to_ascii(filename)\n')
                return
            end
            validateattributes(fname, {'char'}, {'row'}, '', 'fname', 1)
            
            % Open file
            fid = fopen(fname,'wt');
            cleanup = onCleanup(@()fclose(fid));
            
            % Write the header
            fprintf(fid,'# Rotating fluid planet modeled by Concentric Maclaurin Spheroids.\n');
            fprintf(fid,'#\n');
            fprintf(fid,'# Model name: %s\n', obj.name);
            fprintf(fid,'#\n');
            fprintf(fid,'# Scalar quantities:\n');
            fprintf(fid,'# N layers = %d\n',obj.N);
            fprintf(fid,'# Mass M = %g kg\n', double(obj.M));
            fprintf(fid,'# Mean radius       s0 = %0.6e m\n', double(obj.s0));
            fprintf(fid,'# Equatorial radius a0 = %0.6e m\n', double(obj.a0));
            fprintf(fid,'# Rotation parameter m = %0.6f\n', double(obj.mrot));
            fprintf(fid,'# Rotation parameter q = %0.6f\n', double(obj.qrot));
            fprintf(fid,'#\n');
            fprintf(fid,'# Calculated gravity zonal harmonics (x 10^6):\n');
            fprintf(fid,'# J2  = %12.6f\n', obj.J2*1e6);
            fprintf(fid,'# J4  = %12.6f\n', obj.J4*1e6);
            fprintf(fid,'# J6  = %12.6f\n', obj.J6*1e6);
            fprintf(fid,'# J8  = %12.6f\n', obj.J8*1e6);
            fprintf(fid,'#\n');
            fprintf(fid,'# Column data description (MKS):\n');
            fprintf(fid,'# i     - level surface index (increasing with depth)\n');
            fprintf(fid,'# s_i   - mean radius of level surface i\n');
            fprintf(fid,'# a_i   - equatorial radius of level surface i\n');
            fprintf(fid,'# rho_i - density on level surfaces i\n');
            fprintf(fid,'# P_i   - pressure on level surface i\n');
            fprintf(fid,'# m_i   - mass below level surface i\n');
            fprintf(fid,'#\n');
            
            % Write the data
            fprintf(fid,'# Column data:\n');
            fprintf(fid,'# %-4s  ','i');
            fprintf(fid,'%-10s  ','s_i','a_i');
            fprintf(fid,'%-7s  ','rho_i');
            fprintf(fid,'%-10s  ','P_i','m_i');
            fprintf(fid,'\n');
            for k=1:obj.N
                fprintf(fid,'  %-4d  ',k);
                fprintf(fid,'%10.4e  ', double(obj.si(k)));
                fprintf(fid,'%10.4e  ', double(obj.ai(k)));
                fprintf(fid,'%7.1f  ', double(obj.rhoi(k)));
                fprintf(fid,'%10.4e  ', double(obj.Pi(k)), double(obj.mi(k)));
                fprintf(fid,'\n');
            end
        end
    end % End of reporters/exporters block

    %% Private methods
    methods (Access = private)
        function Upu = calc_equipotential_Upu(obj)
            % See ./notes/cms.pdf eqs. 49 and 51
            Upu = NaN(obj.N,1);
            lam = obj.CMS.lambdas;
            til = obj.CMS.JLike.tilde;
            tilp = obj.CMS.JLike.tildeprime;
            tilpp = obj.CMS.JLike.tildeprimeprime;
            q = obj.qrot;
            zet = obj.CMS.zetas(:,2); % 2 for no reason, could pick any angle
            P2k = obj.CMS.Ps.Pnmu(:,2);
            
            n = (0:obj.opts.kmax)';
            xind = obj.CMS.xind;
            xlam = lam(xind);
            nx = length(xind);
            for j=1:obj.N
                U = 0;
                xj = find(xlam <= lam(j), 1);
                if isempty(xj), xj = nx + 1; end % if no xplicit layer below jth
                for i=xj:nx
                    U = U + sum((xlam(i)/lam(j)).^n.*til(i,:)'.*zet(j).^-n.*P2k(:));
                end
                for i=1:xj-1
                    U = U + sum((lam(j)/xlam(i)).^(n+1).*tilp(i,:)'.*zet(j).^(n+1).*P2k(:));
                    U = U + (lam(j)/xlam(i))^3*tilpp(i)*zet(j)^3;
                end
                U = U*(-1/(lam(j)*zet(j)));
                U = U + (1/3)*q*lam(j)^2*zet(j)^2*(1 - P2k(3)); % add rotation
                Upu(j) = U;
            end
        end
    end
    
    %% Access and pseudo-access methods
    methods
        function set.name(obj,val)
            if ~isempty(val)
                validateattributes(val, {'char'}, {'row'})
            end
            obj.name = val;
        end
        
        function set.ai(obj, val)
            assert(isnumeric(val) && isvector(val) && ~any(val<0),...
                'obj.si must be a nonnegative vector.')
            assert(all(diff(val)<=0),'obj.si must be non-ascending.')
            obj.ai = val(:);
        end
        
        function set.rhoi(obj, val)
            assert(isnumeric(val) && isvector(val),...
                'obj.rhoi must be a nonnegative vector.')
            if any(val<0)
                warning('TOFPLANET:assertion','negative density. Is this on purpose?')
            end
            obj.rhoi = val(:);
        end
        
        function val = Pi(obj, ind)
            try
                U = obj.Ui;
                rho = obj.rhoi;
                val = zeros(obj.N,1);
                val(1) = obj.P0;
                val(2:end) = val(1) + cumsum(rho(1:end-1).*diff(U));
                if nargin > 1
                    val = val(ind);
                end
            catch
                val = [];
            end
        end
        
        function val = Ui(obj, ind)
            try
                val = (obj.G*obj.mass/obj.radius)*obj.calc_equipotential_Upu();
                if nargin > 1
                    val = val(ind);
                end
            catch
                val = [];
            end
        end
        
        function set.eos(obj,val)
            if isempty(val)
                obj.eos = [];
                return
            end
            if ~isa(val,'barotropes.Barotrope')
                error('eos must be a valid instance of class Barotrope')
            end
            obj.eos = val(:);
        end
        
        function set.bgeos(obj,val)
            if isempty(val)
                obj.bgeos = [];
                return
            end
            if ~isa(val,'barotropes.Barotrope')
                error('bgeos must be a valid instance of class Barotrope')
            end
            obj.bgeos = val;
        end

        function set.fgeos(obj,val)
            if isempty(val)
                obj.fgeos = [];
                return
            end
            if ~isa(val,'barotropes.Barotrope')
                error('fgeos must be a valid instance of class Barotrope')
            end
            obj.fgeos = val;
        end
        
        function val = get.M(obj)
            try
                drho = [obj.rhoi(1); diff(obj.rhoi)];
                if isempty(obj.si)
                    val = (4*pi/3)*sum(drho.*obj.ai.^3);
                else
                    val = (4*pi/3)*sum(drho.*obj.si.^3);
                end
            catch
                val = [];
            end
        end
        
        function set.mass(obj,val)
            validateattributes(val,{'numeric'},{'positive','scalar'})
            obj.mass = val;
        end
        
        function set.radius(obj,val)
            validateattributes(val,{'numeric'},{'positive','scalar'})
            obj.radius = val;
        end
        
        function val = get.mi(obj)
            % mass _below_ level i
            if isempty(obj.ai) || isempty(obj.rhoi)
                val = [];
            else
                rho = obj.rhoi;
                s = obj.si;
                n = obj.N;
                val(n) = 4*pi/3*rho(n)*s(n)^3;
                for k=n-1:-1:1
                    val(k) = val(k+1) + 4*pi/3*rho(k)*(s(k)^3 - s(k+1)^3);
                end
                val = val';
            end
        end
        
        function val = mzi(obj)
            % heavy element mass below level i
            z = obj.zi;
            if isempty(obj.ai) || isempty(obj.rhoi) || isempty(z)
                val = [];
            else
                val = NaN; %TODO: implement
            end
        end
        
        function val = zi(obj)
            % heavy element mass fraction in layer i
            P = obj.Pi;
            if isempty(obj.bgeos) || isempty(obj.fgeos) || isempty(P)
                val = [];
            else
                roxy = obj.bgeos.density(P);
                roz = obj.fgeos.density(P);
                ro = obj.rhoi;
                val = (1./ro - 1./roxy)./(1./roz - 1./roxy);
                val(~isfinite(val)) = 0;
            end
        end
        
        function val = M_Z(obj)
            try
                val = obj.mzi(1);
            catch
                val = [];
            end
        end
        
        function val = get.N(obj)
            if isempty(obj.ai) || isempty(obj.rhoi)
                val = 0;
            elseif length(obj.ai) == length(obj.rhoi)
                val = length(obj.ai);
            else
                error('length(ai) = %g ~= length(rhoi) = %g',...
                    length(obj.ai),length(obj.rhoi))
            end
        end
        
        function val = get.rhobar(obj)
            if isempty(obj.M) || isempty(obj.s0)
                val = [];
            else
                val = obj.M/(4*pi/3*obj.s0^3);
            end
        end
        
        function val = get.a0(obj)
            if isempty(obj.ai)
                val = [];
            else
                val = obj.ai(1);
            end
        end
        
        function val = get.J2(obj)
            if isempty(obj.Js)
                val = 0;
            else
                val = obj.Js(2);
            end
        end
        
        function val = get.J4(obj)
            if isempty(obj.Js)
                val = 0;
            else
                val = obj.Js(3);
            end
        end
        
        function val = get.J6(obj)
            if isempty(obj.Js)
                val = 0;
            else
                val = obj.Js(4);
            end
        end
        
        function val = get.J8(obj)
            if isempty(obj.Js)
                val = 0;
            else
                val = obj.Js(5);
            end
        end
        
        function val = get.si(obj)
            try
                Vs = NaN(size(obj.CMS.lambdas));
                for j=1:obj.N
                    Vs(j) = obj.CMS.lambdas(j)^3*(obj.CMS.zetas(j,:).^3)*(obj.CMS.gws');
                end
                val = obj.a0*Vs.^(1/3);
            catch
                val = [];
            end
        end
        
        function val = get.s0(obj)
            try
                val = obj.si(1);
            catch
                val = [];
            end
        end
        
        function val = get.mrot(obj)
            try
                val = obj.qrot*(obj.s0/obj.a0)^3;
            catch
                val = [];
            end
        end
        
        function val = NMoI(obj, reduce)
            % C/Ma^2, see eq. 5 in Hubbard & Militzer 2016
            
            if nargin < 2 || isempty(reduce), reduce = 'sum'; end
            reduce = validatestring(reduce, {'sum', 'csum', 'none'});
            
            if isempty(obj.CMS)
                deltas = [obj.rhoi(1); diff(obj.rhoi)];
                lambdas = obj.ai/obj.a0;
                zetas = ones(obj.N, obj.opts.nangles);
            else
                deltas = obj.CMS.deltas;
                lambdas = obj.CMS.lambdas;
                zetas = obj.CMS.zetas;
            end
            [mus, gws] = gauleg(0, 1, obj.opts.nangles); % Abscissas and weights for Gauss-quad
            p2term = 1 - Pn(2, mus);
            num(obj.N) = 0;
            den = 0;
            for j=1:obj.N
                fun1 = deltas(j)*((zetas(j,:)*lambdas(j)).^5).*p2term;
                fun2 = deltas(j)*((zetas(j,:)*lambdas(j)).^3);
                num(j) = fun1*gws';
                den = den + fun2*gws';
            end
            if isequal(reduce, 'none'), val = (2/5)*(num/den); end
            if isequal(reduce, 'sum'), val = (2/5)*sum(num)/den; end
            if isequal(reduce, 'csum'), val = (2/5)*cumsum(num)/den; end
        end
    end % End of access methods block
    
    %% Static methods
    methods (Static)
    end % End of static methods block
end % End of classdef

%% Class-related functions
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
