% Base class for a CONSTANT_PROVIDER object 
% The CONSTANT_PROVIDER will handle getting constants from a source. 
% The source could be an excel file, ini-files, json files or anything else.
% The base class must be sub-classed and extended for each type of source.

classdef CONSTANT_PROVIDER_base_class   
    properties
        source_type % identifier for what type of source is used
        filepath    % path to source file(s)
        const_data  % input from const source
    end
    
    methods
        
        function self = CONSTANT_PROVIDER_base_class(mypath, const_file)
            % CONSTRUCTOR for CONSTANT_PROVIDER_base_class
            %   Must setup and initialize all functionality of the instance.
            %
            %   ARGUMENTS:
            %   mypath:       path to constant file
            %   const_file:   filename of const_file
            %
            %   Must be overloaded when sub-classing.
            
            self.source_type = [];
            self.filepath = [];
        end
        
        function value = get_constant(self, name)
            % GET_CONSTANT  retrieves the value of a particular constant
            %   from the constant source
            %
            %   ARGUMENTS:
            %   name:  the name (string) of the constant to retrieve
            %
            %   RETURNS: 
            %   val:   the the value requested
            
            error('get_constant method not implemented for base_class. You should call a specific sub-class.');
        end
        
        function structure = populate_struct(self, structure)
            % POPULATE_STRUCT  Populates the fields of the provided structure with
            %   values from the constant source
            %
            %   ARGUMENTS:
            %   structure:  a structure with empty fields to be populated
            %
            %   RETURNS: 
            %   structure:  the input structure with fields populated.
            %
            %   Must be overloaded when sub-classing.
            
            error('populate struct method not implemented for base_class. You should call a specific sub-class.');
        end

    end
end

