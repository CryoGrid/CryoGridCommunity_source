
classdef PARAMETER_PROVIDER_YAML < PARAMETER_PROVIDER_base_class
    % Child class of PARAMETER_PROVIDER
    % The PARAMETER_PROVIDER_YAML will handle getting model
    % parameters from a YAML file.
    
    properties
        config_file
    end
    
    
    methods
        
        function self = PARAMETER_PROVIDER_YAML(mypath, config_file)
            % CONSTRUCTOR for PARAMETER_PROVIDER_YAML
            %   Reads in parameters from config file of type 'yml'.
            %
            %   ARGUMENTS:
            %   mypath:       path to config file
            %   config_file:  filename of config_file

            self.config_file = config_file;
            self.config_data = ReadYaml([mypath '/' config_file]);
            self.source_type = 'yml';
            self.filepath = mypath;
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

            id = [];
            
            for i = 1:length(self.config_data.(section))
                if strcmp(self.config_data.(section){1,i}.type, name) & ...
                        self.config_data.(section){1,i}.index == index
                    id = i;
                end
            end
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
            
            name = self.config_data.(section){1,id}.type;
            index = self.config_data.(section){1,id}.index;
        end
       
        function structure = populate_struct(self, structure, section, name, index)
            % POPULATE_STRUCT  Populates the fields of the provided structure with
            %   values from the parameter source (here xlsx file)
            %
            %   NOTICE: The index is here the actual index specified for a
            %     certain class in the xlsx file.
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
            
            % Get list of fieldnames in the requested structure
            fn = fieldnames(structure);
            id = self.get_class_id_by_name_and_index(section, name, index);
            
            assigned_fields = [];
            
            for i = 1:size(fn,1)
                if any(strcmp(fn{i}, fieldnames(self.config_data.(section){1,id})))
                    if isstruct(self.config_data.(section){1,id}.(fn{i}))
                        fn_substruct = fieldnames(structure.(fn{i}));
                        assigned_subfields = [];
                        assigned_fields = [assigned_fields fn(i)];
                        for j = 1:size(fn_substruct,1)
                            if any(strcmp(fn_substruct{j},fieldnames(self.config_data.(section){1,id}.(fn{i}))))
                                structure.(fn{i}).(fn_substruct{j}) = self.config_data.(section){1,id}.(fn{i}).(fn_substruct{j})';
                                assigned_subfields = [assigned_subfields fn_substruct(j)];                       
                            end
                        end
                    else
                        structure.(fn{i}) = self.config_data.(section){1,id}.(fn{i});
                        assigned_fields = [assigned_fields fn(i)];
                    end
                end
            end
            
            for j = 1:size(fn,1)
                    if ~any(strcmp(assigned_fields, fn{j}))
                    %     if one of the structure field is not populated, displays a warning message that "fieldname"
                    %     was not populated
                        warning(['The structure field ' fn{j} ' was not populated']);
                    end
            end
            
            if exist ('fn_substruct')
                for k = 1:size(fn_substruct,1)
                    if ~any(strcmp(assigned_subfields, fn_substruct{k}))
                    %   if not populated, give a warning message that
                    %       "fieldname.substructure fieldname" was not populated
                        warning(['The sub-structure field ' fn_substruct{k} ' was not populated']);
                    end
                end
            end
        end
        
    end
end




