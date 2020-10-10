
classdef PARAMETER_PROVIDER_EXCEL < PARAMETER_PROVIDER_base_class
    % Child class of PARAMETER_PROVIDER
    % The PARAMETER_PROVIDER_EXCEL will handle getting model
    % parameters from an Excel file.
    
    properties
        config_file
    end
    
    
    methods
        
        function self = PARAMETER_PROVIDER_EXCEL(mypath, config_file)
            % CONSTRUCTOR for PARAMETER_PROVIDER_EXCEL
            %   Reads in parameters  from config file of type 'xlsx'.
            %
            %   ARGUMENTS:
            %   mypath:       path to config file
            %   config_file:  filename of config_file

            self.config_file = config_file;
            self.config_data = read_excel2cell([mypath '/' config_file]);
            self.source_type = 'xlsx';
            self.filepath = mypath;
            self.tile_info = self.get_tile_information('TILE_IDENTIFICATION');
        end
		
        function forcing_file = get_forcing_file_name(self, section) 
		% GET_FORCING_FILE_NAME  Get forcing file name from specified section of 
			% the configuration file.
            %
            %   ARGUMENTS:
            %   section: a string specifying the forcing section to list 
            %                (e.g. 'FORCING')
            %
            %   RETURNS: 
            %   forcing_file: string with forcing file name extracted from the 
			%   				configuration file
			
            pos_list = self.get_range_section(section);

            for i=1:size(pos_list,1)
                % get cells defining the current class
                section_data = self.config_data(pos_list(i,1):pos_list(i,2),:);
                if ~isempty(find(strcmp(section_data(:,1:2),'filename')))
					row_index = find(strcmp(section_data(:,1:2),'filename'));
                    if ~isnan(section_data{row_index,2})
						forcing_file = section_data{row_index,2};
					else
						error('The name of the forcing file is not specified in the configuration file')
					end 
				end
			end
		end

		function [tile_info] = get_tile_information(self, section) 
		% GET_TILE_INFORMATION Get tile information and inputs for the tile 
			% builder from specified section of the configuration file.
            %
            %   ARGUMENTS:
            %   section: a string specifying the tile identification section  
            %                (e.g. 'TILE_IDENTIFICATION')
            %
            %   RETURNS: 
            %   tile_info: a structure defining the tile type, location, class 
			%		combination, etc.
			
			tile_info = struct();
            pos_list = self.get_range_section(section);
            section_data = self.config_data(pos_list(1,1):pos_list(1,2),:);
			
            for i=2:size(section_data,1)-1
                if strcmp(section_data{i,2}, 'LIST') 
                    list_end = 0;
                    for j=3:size(section_data,2)
                        if strcmp(section_data{i,j}, 'END') 
                            list_end = j-1;
                            break;
                        end
                    end
                    if list_end > 0
                        tile_info.(section_data{i,1}) = section_data(i,3:list_end);       
                    else
                        warning(['End of list not identified for field ' section_data{i,1} ' in TILE_IDENTIFICATION section']);
                    end
                else
                    tile_info.(section_data{i,1}) = section_data{i,2};
                end
		    end	
			
			tile_info.coordinates = str2num(tile_info.coordinates);
			
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
            %                class_name, class_index, from_row, to_row
            %                (last two columns indicate location in cell 
            %                 array read from xlsx file)
            
            % get range of cells in cell array that make up the section
            % and the individual classes defined in the section.
            
            pos_list = self.get_range_section(section);
            %pos_list = self.get_range(section);
            
            % define the table to populate
%             class_list = table('size', [size(pos_list,1),4], ...
%                               'VariableTypes', {'string', 'uint32', 'uint8', 'uint8'}, ...
%                               'VariableNames',
%                               {'class_name','class_index','from_row','to_row'});
%                               %chnged from uiont8 to uint32, S.Westermann, Oct 2020
            class_list = table('size', [size(pos_list,1),4], ...
                              'VariableTypes', {'string', 'uint32', 'uint32', 'uint32'}, ...
                              'VariableNames', {'class_name','class_index','from_row','to_row'});
                          
            % loop over all class definitions 
            % (Adapted to cases where section_data does not include the 
            % header 'index' in column 2. Happen when section keywords like 
            % STRAT_linear, STRAT_classes and STRAT_layers are used, and the headers are dismissed). 
            for i=1:size(pos_list,1)
                % get cells defining the current class
                section_data = self.config_data(pos_list(i,1):pos_list(i,2),:);
                if ~any(strcmp(section_data(:,1:2),'index'))
                    class_list.class_name(i) = section_data{1,1};
                    class_list.class_index(i) = section_data{1,2};
                else
                    % extract name, index and from and to cell ids
                    class_list.class_name(i) = section_data{2,1};
                    class_list.class_index(i) = section_data{2,2};
                end
                class_list.from_row(i) = pos_list(i,1);
                class_list.to_row(i) = pos_list(i,2);
            end
            
            if size(class_list, 1) == 0
                error(['The section "' section '" was not found in the configuration file!']);
            end
            
        end
        
        
        function id = get_class_id_by_name_and_index(self, section, name, index)
            % GET_CLASS_ID_BY_NAME_AND_INDEX  Gets the class id of a class
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
                        
            class_list = self.get_class_list(section);
            id = find(strcmp(class_list.class_name, name) & ...
                             class_list.class_index == index);
            if isempty(id)
                id = [];  % ensure always same type if not found
            end
            
        end
        
        
        function [name, index] = get_class_name_and_index_by_id(self, section, id)
            % GET_CLASS_NAME_AND_INDEX_BY_ID  Gets the class_name of the id'th class
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
            
            class_list = self.get_class_list(section);
            name = class_list.class_name(id);
            index = class_list.class_index(id);
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
            
            if ~isstruct(structure)
                % Do nothing if we are not passed a proper structure.
                % Happens e.g. when PARA has no fields (empty)
                return
            end
            
            % Get list of fieldnames in the requested structure
            fn = fieldnames(structure);
            % Extract relevant section from Excel file data
            cl = self.get_class_list(section);
            from_row = cl.from_row(cl.class_name==name & cl.class_index==index);
            to_row = cl.to_row(cl.class_name==name & cl.class_index==index);
            section_data = self.config_data(from_row:to_row,:);
            
            assigned_fields = [];
            
            for i = 1:size(section_data,1)
            % Loop over section_data
                if strcmp(section_data{i,1},'TOP')
                % If TOP / BOTTOM is found:
                %    Extract table section from row TOP-2 to line BOTTOM
                %    Remove these lines from section_data
                %    Store in TABLE_DATA temporary matlab table
                %        using the provided column names in row TOP-2
            
                    pos_list = get_range_TOP_BOTTOM(section_data);
                
                    %variable_names = {};
                    %variable_values = [];
                    table_data = table();
                    x = 1;
                    field = cell2mat(section_data(pos_list(1,1)-3,x));
                    while x <= size(section_data,2) && ~isnan(field(1))
                        %field = cell2mat(section_data(pos_list(1,1)-3, x)) %CHANGED Sebastian Oct 2020!
                        
                        table_data.(section_data{pos_list(1,1)-3, x}) = section_data(pos_list(1,1):pos_list(1,2), x);
                        x = x+1;
                        %field = cell2mat(section_data(pos_list(1,1)-3, x));  
                        
                         if x <= size(section_data,2)                       %CHANGED Sebastian Oct 2020!
                             field = cell2mat(section_data(pos_list(1,1)-3, x));
                         end
                    end
                    for j = 1:size(fn,1) 
                    % Loop over all fieldnames in fn
                        if isstruct(structure.(fn{j}))
                        %   if fieldname is a substructure:
                        %     Loop over all fields in substructure
                        %       find substructure fieldname column in TABLE_DATA
                        %       Assign this column data to structure.substructure.fieldname
                            assigned_fields = [assigned_fields fn(j)];
                            fn_substruct = fieldnames(structure.(fn{j}));
                            assigned_subfields = [];
							if ~isempty(fn_substruct)
                                for k = 1:size(fn_substruct,1)
									if any(strcmp(table_data.Properties.VariableNames, fn_substruct{k}))
										structure.(fn{j}).(fn_substruct{k}) = table_data.(fn_substruct{k});
										assigned_subfields = [assigned_subfields fn_substruct(k)];
									end
								end
							else
							    for k = 1:length(table_data.Properties.VariableNames)
										structure.(fn{j}).(table_data.Properties.VariableNames{k}) = table_data.(table_data.Properties.VariableNames{k});
								end	
							end
                        end
                    end
                elseif any(strcmp(fn, section_data{i,1}))
                %     find "fieldname" parameter in section data and assign to
                %       structure.fieldname
                    structure.(section_data{i,1}) = section_data{i,2};
                    assigned_fields = [assigned_fields section_data(i,1)];
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


