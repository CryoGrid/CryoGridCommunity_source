% Base class for a PARAMETER_PROVIDER object 
% The PARAMETER_PROVIDER will handle getting model parameters from a source. 
% The source could be an excel file, ini-files, json files or anything else.
% The base class must be sub-classed and extended for each type of source.

classdef PARAMETER_PROVIDER_base_class
    properties
        source_type % identifier for what type of source is used
        filepath    % path to source file(s)
        config_data % input from config source
    end
    
    
    methods
        
        %mandatory functions for each class

        function self = PARAMETER_PROVIDER_base_class(mypath, config_file)
            % CONSTRUCTOR for PARAMETER_PROVIDER_base_class
            %   Must setup and initialize all functionality of the instance.
            %
            %   ARGUMENTS:
            %   mypath:       path to config file
            %   config_file:  filename of config_file
            %
            %   Must be overloaded when sub-classing.
            
            self.source_type = [];
            self.filepath = [];
        end

        
        function class_list = get_class_list(self, section)
            % GET_CLASS_LIST  Gets a table of available classes in the
            %   specified section.
            %
            %   ARGUMENTS:
            %   section: a string specifying which section to list 
            %                (e.g. 'FORCING' or 'OUT')
            %
            %   RETURNS: 
            %   class_list: a table with the following columns: 
            %                class_name, class_index
            %                (more columns may be added if needed by sub-class)
            %
            %   Must be overloaded when sub-classing.
           
            error('get_class_list method not implemented for base_class.');
        end

        
        function id = get_class_id_by_name_and_index(self, section, name, index)
            % GET_CLASS_ID_BY_NAME_AND_INDEX  Get the class id of a class
            %   given by its class name and class index.
            %
            %   NOTICE: The id here is not the class index, but the order 
            %      in which it is defined in the excel file.
            %      The index is the index specified in the source in the
            %      corresponding class definition.
            %
            %   ARGUMENTS:
            %   section:  a string specifying which section to list 
            %                (e.g. 'FORCING' or 'OUT')
            %   name:     name of the class to find
            %   index:    defined index of the class to find
            %
            %   RETURNS: 
            %   id:       the id of the class (the position in the section) 
            %
            %   Must be overloaded when sub-classing.
            
            error('get_class_name_by_id method not implemented for base_class.');
        end
          
        
        function [name, index] = get_class_name_and_index_by_id(self, section, id)
            % GET_CLASS_NAME_AND_INDEX_BY_ID  Get the class_name of the id'th class
            %   defined in the given section.
            %
            %   NOTICE: The id here is not the class index, but the order 
            %      in which it is defined in the excel file.
            %      If we ask for id=1, we will get e.g. the first OUT 
            %      class defined in the OUT section 
            %      (since more OUT classes may be defined)
            %
            %   ARGUMENTS:
            %   section:  a string specifying which section to list 
            %                (e.g. 'FORCING' or 'OUT')
            %   id:       the id of the class (the position in the section) 
            %
            %   RETURNS: 
            %   class_name:   string
            %   class_index:  integer
            %
            %   Must be overloaded when sub-classing.
            
            error('get_class_name_and_index_by_id method not implemented for base_class.');
        end
        
        
        function structure = populate_struct(self, structure, section, name, index)
            % POPULATE_STRUCT  Populates the fields of the provided structure with
            %   values from the parameter source
            %
            %   ARGUMENTS:
            %   structure:  a structure with empty fields to be populated
            %   section: a string specifying which section to adress
            %                (e.g. 'FORCING' or 'OUT')
            %   name:    name of the class for which to obtain parameters
            %   index:   index of class for which to obtain parameters
            %
            %   RETURNS: 
            %   structure:  the input structure with fields populated.
            %
            %   Must be overloaded when sub-classing.
            
            error('get_parameters method not implemented for base_class. You should call a specific sub-class.');
        end
    end
end

