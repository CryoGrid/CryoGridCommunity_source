%========================================================================
% CryoGrid BASE class for all LATERAL classes 
% contains variables and initialization routines
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, Oct 2020
%========================================================================


classdef BASE_LATERAL < matlab.mixin.Copyable
    
    properties
        class_index
        CONST    %constants
        PARA     %external service parameters, all other
        STATVAR  %energy, water content, etc.
        TEMP     %derivatives in prognostic timestep and optimal timestep
        PARENT
    end
    
    methods
        
%         function self = BASE_LATERAL(index, pprovider, cprovider)
%             % CONSTRUCTOR for BASE_LATERAL class
%             %   Reads in class parameters from parameter providers
%             %
%             %   ARGUMENTS:
%             %   index:      user defined class index
%             %   pprovider:  instance of PARAMETER_PROVIDER class
%             %   cprovider:  instance of CONSTANT_PROVIDER class
%             
%             if index ~=-1
%                 self.class_index = index;
%                 self = self.initialize();
%                 self = self.populate_PARA(pprovider);
%                 self = self.populate_CONST(cprovider);
%             end
%         end

        function self = initialize_excel(self)
            
        end        


        function self = initialize(self)
            % INITIALIZE  Initializes all properties needed by the class.
            self = self.initialize_PARA();
            self = self.initialize_TEMP();
            self = self.initialize_STATVAR();
            self = self.initialize_CONST();
        end
        
        function self = initialize_PARA(self)
            % INITIALIZE_PARA  Initializes PARA structure. Only the parameters populated by the PARAMETER_PROVIDER should be defined here, additional parameters will be dinamically created during the final stage of the initialization (finalize_setup).
            self = self.provide_PARA();
        end
        
        function self = initialize_TEMP(self)
            % INITIALIZE_TEMP  Initializes TEMP structure.
            self.TEMP = struct();
        end
        
        function self = initialize_STATVAR(self)
            % INITIALIZE_STATVAR  Initializes STATVAR structure.
            self = self.provide_STATVAR();
        end
        
        function self = initialize_CONST(self)
            % INITIALIZE_CONST  Initializes CONST structure. Only the parameters populated by the CONSTANT_PROVIDER should be defined here.
            self = self.provide_CONST();
        end
        
        function self = populate_PARA(self, pprovider)
            % POPULATE_PARA  Updates the PARA structure with values from pprovider.
            %
            %   ARGUMENTS:
            %   pprovider:  instance of PARAMETER_PROVIDER class
            
            self.PARA = pprovider.populate_struct(self.PARA, 'LATERAL_IA_CLASS', class(self), self.class_index);
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

