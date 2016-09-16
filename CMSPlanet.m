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
        cms     % a CMS object
        eos     % a barotrope object
        opts    % struct with opts, both CMS and CMSPlanet
    end
    properties (Dependent)
        qrot    % rotation parameter 
        nlayers % layers of constant density 
        Jn      % external gravity moments
        s0      % mean radius
        b0      % polar radius
        f0      % flattening, a.k.a, oblateness: (a - b)/a
    end
    properties (Access = private)
        si
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
            
            % Save opts struct and call InitPlanet to finish the construction
            obj.opts = op; % (calls set.opts and may have side effects!)
            obj.InitPlanet(op);
        end
    end % End of constructor block
    
    %% Public methods
    methods (Access = public)
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
        function set.opts(obj,val)
            obj.opts = val;
            obj.cms.opts = val; %#ok<MCSUP>
        end
        
        function val = get.nlayers(obj)
            val = obj.cms.nlayers;
        end
        
        function set.nlayers(obj,val)
            obj.cms.nlayers = val;
        end
        
        function val = get.qrot(obj)
            val = obj.cms.qrot;
        end
        
        function set.qrot(obj,val)
            obj.cms.qrot = val;
        end
        
        function val = get.Jn(obj)
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
    end % End of access methods block
    
    %% Static methods
    methods (Static)
        
    end % End of static methods block
end % End of classdef

