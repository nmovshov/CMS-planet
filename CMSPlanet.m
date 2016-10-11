classdef CMSPlanet < handle
    %CMSPLANET Interior model of rotating fluid planet.
    %   This class implements a model of a rotating fluid planet using the
    %   Concentric Maclaurin Spheroids technique. A CMSPlanet object is defined
    %   by a given mass, equatorial radius, rotation period, and barotrope.
    
    %% Properties
    properties (Access = public)
        name    % model name
        desc    % model description (one line)
        a0      % equatorial radius
        M       % total mass
    end
    properties (Dependent)
        ai      % layer equatorial radii
        bi      % layer polar radii
        si      % layer mean radii
        rhoi    % layer densities
        Pi      % layer surface pressure (top of layer)
    end
    properties
        cms     % a CMS object
        eos     % a barotrope object
    end
    properties (Dependent)
        opts    % struct with opts, both CMS and CMSPlanet
        qrot    % rotation parameter 
        nlayers % layers of constant density 
        Js      % external gravity moments (even, 0:2:kmax degree)
        s0      % mean radius
        b0      % polar radius
        f0      % flattening, a.k.a, oblateness: (a - b)/a
        rho0    % reference density (uses equatorial radius)
        rho_s   % mean density (uses mean radius)
        M_calc  % mass from current state of cms
    end
    properties (Access = private)
        
    end
    
    %% The constructor
    methods
        function obj = CMSPlanet(varargin)
            %CMSPLANET Class constructor.
            
            % Empty, integer, or name/value arguments.
            warning off CMS:obsolete
            if nargin == 0
                op = cmsset();
            else
                if isnumeric(varargin{1})
                    varargin = {'nlayers',varargin{1},varargin{2:end}};
                end
                op = cmsset(varargin{:});
            end
            warning on CMS:obsolete
            
            % Call InitPlanet to finish the construction
            obj.InitPlanet(op);
        end
    end % End of constructor block
    
    %% Public methods
    methods (Access = public)
        function ET = relax_to_HE(obj)
            ET = obj.cms.relax();
        end
        
        function fac = match_total_mass(obj)
            % Rescale layer densities to match planet mass; return scale factor.
            if isempty(obj.rhoi), fac = []; return, end
            fac = obj.M/obj.M_calc;
            obj.rhoi = obj.rhoi*fac;
        end
        
        function ah = show(obj)
            % Visualize a CMSPlanet object, return axes handle.
            
            % Require R2016a to use the amazing polarplot features
            if verLessThan('matlab','9')
                error('CMS plotting requires R2016a or later')
            end
            
            ah = obj.cms.plot();
        end
    end % End of public methods block
    
    %% Private methods
    methods (Access = private)
        function InitPlanet(obj,op)
           % (Re)Initialize a CMSPlanet object.
           
           obj.cms = ConcentricMaclaurinSpheroids(op);
        end
    end % End of private methods block
    
    %% Access methods
    methods
        function val = get.opts(obj)
            val = obj.cms.opts;
        end
        
        function set.opts(obj,val)
            obj.cms.opts = val;
        end
        
        function val = get.nlayers(obj)
            val = obj.cms.nlayers;
        end
        
        function set.nlayers(~,~)
            msg = ['Changing number of layers in an existing CMSPlanet ',...
                'makes no sense; create a new object instead.'];
            error(msg)
        end
        
        function val = get.qrot(obj)
            val = obj.cms.qrot;
        end
        
        function set.qrot(obj,val)
            try
                (val == 1); %#ok<EQEFF>
                obj.cms.qrot = double(val);
            catch ME
                if (strfind(ME.message,'Attempt to compare')) == 1
                    error('Check dimensions of qrot')
                else
                    rethrow(ME)
                end
            end
        end
        
        function val = get.Js(obj)
            val = obj.cms.Jn;
        end
        
        function set.M(obj,val)
            try
                assert(isnumeric(val))
                assert(isscalar(val))
                assert(double(val) > 0)
                obj.M = val;
            catch
                error('Mass must be a positive scalar.')
            end
        end
        
        function set.a0(obj,val)
            try
                assert(isnumeric(val))
                assert(isscalar(val))
                assert(double(val) > 0)
                obj.a0 = val;
            catch
                error('Radius must be a positive scalar.')
            end
        end
        
        function val = get.ai(obj)
            val = [];
            if ~isempty(obj.a0)
                val = obj.a0*obj.cms.lambdas;
            end
        end
        
        function set.ai(obj,val)
            validateattributes(val,{'numeric'},{'real','finite','vector'},...
                '','ai')
            assert(numel(val) == obj.nlayers,...
                'length(ai) = %g ~= nlayers = %g',...
                numel(val),obj.nlayers)
            assert(all(double(val) > 0), 'layer radii must be positive')
            obj.a0 = val(1);
            obj.cms.lambdas = val/obj.a0;
        end
        
        function val = get.rhoi(obj)
            val = [];
            if ~isempty(obj.rho0)
                val = cumsum(obj.cms.deltas);
                val = val*obj.rho0;
            end
        end
        
        function set.rhoi(obj,val)
            validateattributes(val,{'numeric'},{'real','finite','vector'},...
                '','rhoi')
            assert(numel(val) == obj.nlayers,...
                'length(rhoi) = %g ~= nlayers = %g',...
                numel(val),obj.nlayers)
            assert(all(double(val) >= 0), 'layer densities must be nonnegative')
            assert(~isempty(obj.a0) && ~isempty(obj.M),...
                'cannot set layer densities before setting mass and radius')
            obj.cms.deltas = [val(1); diff(val(:))]/obj.rho0;
        end
        
        function set.eos(obj,val)
            if ~isa(val,'barotropes.Barotrope')
                error('eos must be a valid instance of class Barotrope')
            end
            obj.eos = val;
        end
        
        function val = get.s0(obj)
            val = obj.a0*obj.cms.ss(1);
        end
        
        function val = get.b0(obj)
            val = obj.a0*obj.cms.bs(1);
        end
        
        function val = get.f0(obj)
            val = (obj.a0 - obj.b0)/obj.a0;
        end
        
        function val = get.rho0(obj)
            if isempty(obj.M) || isempty(obj.a0)
                val = [];
            else
                val = obj.M./(4*pi/3*obj.a0^3);
            end
        end

        function val = get.rho_s(obj)
            if isempty(obj.M) || isempty(obj.s0)
                val = [];
            else
                val = obj.M./(4*pi/3*obj.s0^3);
            end
        end

        function val = get.M_calc(obj)
            if isempty(obj.rhoi), val = []; return, end
            drho = [obj.rhoi(1); diff(obj.rhoi)];
            val = (4*pi/3)*(obj.a0^3)*sum(drho.*obj.cms.Vs);
        end
        
        function val = get.bi(obj)
            val = [];
            if ~isempty(obj.a0)
                val = obj.a0*obj.cms.bs;
            end
        end
        
        function val = get.si(obj)
            val = [];
            if ~isempty(obj.a0)
                val = obj.a0*obj.cms.ss;
            end
        end
        
        function val = get.Pi(obj)
            if isempty(obj.rho0), val = []; return, end
            try
                si = setUnits; % if you have physunits
            catch
                si = setFUnits; % if you don't have physunits
            end
            G = si.gravity;
            U = mean(obj.cms.Upu, 2)*G*obj.M_calc/obj.a0;
            val = zeros(obj.nlayers, 1)*si.Pa;
            for j=2:obj.nlayers
                val(j) = val(j-1) + obj.rhoi(j-1)*(U(j) - U(j-1));
            end
        end
        
        function set.Pi(~,~)
            error('Pi is a calculated property and cannot be assigned')
        end
    end % End of access methods block
    
    %% Static methods
    methods (Static)
        
    end % End of static methods block
end % End of classdef
