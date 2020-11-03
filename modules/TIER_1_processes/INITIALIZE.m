%========================================================================
% CryoGrid TIER1 library class for functions related to GROUND and LATERAL_IA class initialization
% interfaces with PROVIDER classes
% T. Ingeman-Nilsen, J. Scheer, S. Westermann, October 2020
%========================================================================

classdef INITIALIZE < BASE
    
    
    methods
        function self = INITIALIZE(index, pprovider, cprovider, forcing)               % Temporary definition, to allow old code to run
            %function self = GROUND_base_class(index, pprovider, cprovider, forcing)      % Definition to be used when old code is no longer supported
            % CONSTRUCTOR for GROUND_base_class
            %   Reads in class parameters from the specified file.
            %
            %   ARGUMENTS:
            %   index:      user defined class index
            %   pprovider:  instance of PARAMETER_PROVIDER class
            %   cprovider:  instance of CONSTANT_PROVIDER class
            %   forcing: 	instance of FORCING class
            
            % The following is only needed to allow legacy code to run
            % May be removed when deprecated functions are removed
            %             nargin
            %             if nargin==4
            %                 index = varargin{1};
            %                 pprovider = varargin{2};
            % 				cprovider = varargin{3};
            % 				forcing = varargin{4};
            %             else
            %                 st = dbstack;
            %                 warning(['DEPRECATION WARNING: Instantiating ' st.name '() with no arguments is deprecated.' newline,...
            %                          'You should update your code to take advantage of new IO interface.']);
            %                 return
            %             end
            % End allow legacy code
            if index ~=-1
                self.class_index = index;
                self = self.initialize();
                self = self.populate_PARA(pprovider);
                self = self.populate_CONST(cprovider);
                %self = self.finalize_setup(forcing);
            end
        end
        
        function self = initialize(self)
            % INITIALIZE  Initializes all properties needed by the class.
            
            self.PREVIOUS = [];
            self.NEXT = [];
            self.IA_PREVIOUS = [];
            self.IA_NEXT = [];
            self = self.initialize_PARA();
            self = self.initialize_TEMP();
            self = self.initialize_STATVAR();
            self = self.initialize_CONST();
        end
        
        function self = initialize_PARA(self)
            % INITIALIZE_PARA  Initializes PARA structure. Only the parameters populated by the PARAMETER_PROVIDER should be defined here, additional parameters will be dinamically created during the final stage of the initialization (finalize_setup).
            
            %self = self.provide_PARA();
            self = provide_PARA(self);
        end
        
        function self = initialize_TEMP(self)
            % INITIALIZE_TEMP  Initializes TEMP structure.
            
            self.TEMP = struct();
        end
        
        function self = initialize_STATVAR(self)
            % INITIALIZE_STATVAR  Initializes STATVAR structure.
            
            %self = self.provide_STATVAR();
            self = provide_STATVAR(self);
        end
        
        function self = initialize_CONST(self)
            % INITIALIZE_CONST  Initializes CONST structure. Only the parameters populated by the CONSTANT_PROVIDER should be defined here.
            
            %self = self.provide_CONST();
            self = provide_CONST(self);
        end
        
        function self = populate_PARA(self, pprovider)
            % POPULATE_PARA  Updates the PARA structure with values from pprovider.
            %
            %   ARGUMENTS:
            %   pprovider:  instance of PARAMETER_PROVIDER class
            
            %self.PARA = pprovider.populate_struct(self.PARA, 'GROUND_CLASS', class(self), self.class_index);  %CHANGED_SEBASTIAN!!!
            self.PARA = pprovider.populate_struct(self.PARA, 'SUBSURFACE_CLASS', class(self), self.class_index);  %CHANGED_SEBASTIAN!!!
        end
        
        function self = populate_CONST(self, cprovider)
            % POPULATE_CONST  Updates the CONST structure with values from cprovider.
            %
            %   ARGUMENTS:
            %   cprovider:  instance of CONSTANT_PROVIDER class
            
            self.CONST = cprovider.populate_struct(self.CONST);
        end
        

        
    end
end

