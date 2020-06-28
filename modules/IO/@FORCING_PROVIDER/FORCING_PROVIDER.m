% Base class for a FORCING_PROVIDER object 
% The FORCING_PROVIDER will handle getting forcing data from a source. 
% The source could be a mat file or anything else (!! only handles mat files at this stage !!).
% The base class must be sub-classed and extended for each type of source.

classdef FORCING_PROVIDER
    properties
        file_path
        forcing_file
        DATA
%        const_file % Is it necessary ?
    end
    
    
    methods
        
        function self = FORCING_PROVIDER(mypath, forcing_file)
            % CONSTRUCTOR for FORCING_PROVIDER
            %   Reads in forcing data from the specified file.
            %
            %   ARGUMENTS:
            %   mypath:        path to forcing dataset file
            %   forcing_file:  filename of forcing_file
            
            self.file_path = mypath;
            self.forcing_file = forcing_file;
            self = self.load_forcing_from_mat();
        end


        function self = load_forcing_from_mat(self)
            % LOAD_FORCING_FROM_MAT  Loads the data from source and
            %   and asigns it to DATA property.

            tmp = load(fullfile(self.file_path, self.forcing_file), 'FORCING');
            self.DATA = tmp.FORCING.data;
        end
        
        
        function structure = populate_struct(self, structure)
            % POPULATE_STRUCT  Populates the fields of the provided struct with
            %   data from the forcing source 
            %
            %   ARGUMENTS:
            %   structure:  a structure with empty fields to be populated
            %
            %   RETURNS: 
            %   structure:  the input structure with fields populated.            
            
            % Get list of fieldnames in the requested structure
            fn = fieldnames(structure);
            
            % Initialize cell array to hold names of fields to which values have been assigned
            assigned_fields = cell(size(fn));
            
            % Now loop over all field names in structure
            for i = 1:length(fn)
                if isfield(self.DATA, fn{i})
                    structure.(fn{i}) = self.DATA.(fn{i});
                    
                    % Add field name to assigned_fields
                    assigned_fields{i} = fn{i};
                end                    
            end
            
            % assigned_fiels is currently not used
            
            % It keeps track of all fieldnames in input structure
            % to which values are assigned.
            % It could be used to warn user that some fields were not
            % assigned.
        end
    end
end
            
            