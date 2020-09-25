
classdef CONSTANT_PROVIDER_EXCEL < CONSTANT_PROVIDER_base_class
    % Child class of CONSTANT_PROVIDER
    % The CONSTANT_PROVIDER_EXCEL will handle getting constants from an Excel file.
    
    properties
        const_file
    end
    
    methods
        
        function self = CONSTANT_PROVIDER_EXCEL(mypath, const_file)
            % CONSTRUCTOR for CONSTANT_PROVIDER_EXCEL
            %   Reads in constants from constant file of type 'xlsx'.
            %
            %   ARGUMENTS:
            %   mypath:       path to constant file
            %   const_file:   filename of const_file

            self.const_file = const_file;
            self.const_data = read_excel2cell([mypath '/' const_file]);
            self.source_type = 'xlsx';
            self.filepath = mypath;
        end
        
         
        function structure = populate_struct(self, structure)
            % POPULATE_STRUCT  Populates the fields of the provided structure with
            %   values from the constant source (here xlsx file)
            %
            %   NOTICE: The index is here the actual index specified for a
            %     certain class in the xlsx file.
            %
            %   ARGUMENTS:
            %   structure:  a structure with empty fields to be populated
            %
            %   RETURNS: 
            %   structure:  the input structure with fields populated.
            
            if ~isstruct(structure)
                % Do nothing if we are not passed a proper structure.
                % Happens e.g. when CONST has no fields (empty)
                return
            end
            
            % Get list of fieldnames in the requested structure
            fn = fieldnames(structure);
            
            % Extract relevant section from Excel file data
            section_data = self.const_data; 
            
            % Initialize cell array to hold names of fields to which values have been assigned
            assigned_fields = cell(size(fn));
            
            % Loop over all rows, and extract variables
            for i = 1:size(section_data,1)
                if any(strcmp(fn, section_data{i,1}))
                    % Current cell of first column matches fieldname in structure
                    structure.(section_data{i,1}) = section_data{i,2};                    
                    % Add field name to assigned_fields
                    assigned_fields{i} = section_data{i,1};
                end
            end
            
            % Now loop over all field names in structure
            for i = 1:length(fn)
                if ~any(strcmp(assigned_fields, fn{i}))
                    % if this field was not yet assigned, look for a table that includes the field
                    st = dbstack;
                    warning(['!!! Table lookup in ' st.name ' is not yet implemented.']);
                end
            end                                 
        end
    end
end

