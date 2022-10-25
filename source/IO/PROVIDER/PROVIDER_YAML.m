classdef PROVIDER_YAML < BASE_PROVIDER

    methods
        
        function provider = assign_paths_yaml(provider, run_name, result_path, constant_file)
            
			%warning('WARNING: YAML provider not fully implemented and tested!')
			
            constant_file = [result_path run_name '/' constant_file '.yml'];
            parameter_file = [result_path run_name '/' run_name '.yml'];

            provider.PARA.run_name = run_name;
            provider.PARA.result_path = result_path;
            provider.PARA.parameter_file = parameter_file;
            provider.PARA.constant_file = constant_file;
            %provider.PARA.forcing_path = forcing_path;
        end
       
        
        function provider = read_const_yaml(provider)
            
            data = yaml.ReadYaml(provider.PARA.constant_file);
            fnames = fieldnames(data.CONST);
            for i=1:length(fnames)
                if ~isnan(data.CONST.(fnames{i}))
                    provider.CONST.(fnames{i}) = data.CONST.(fnames{i});
                end
            end
        end
        
        function provider = read_parameters_yaml(provider)
            
            %read the parameter file
            data = yaml.ReadYaml(provider.PARA.parameter_file);

            % get all the CLASS sections defined in the file
            class_sections = fieldnames(data); 
            
            % loop over all sections
            for i = 1:length(class_sections)
                section = data.(class_sections{i});
                
                if isempty(section)
                    % if section is set to 'null' in yaml file
                    % skip processing it
                    continue
                end

                % loop over definitions in section
                for k = 1:length(section)
                    % instantiate new class
                    class_name = section{k}.name;
                    class_index = section{k}.index;
                    class_func = str2func(class_name);
                    new_class = class_func();  %initialize the class
                    
                    %initialize CONST
                    new_class = provide_CONST(new_class);
                    if ~isempty(new_class.CONST)
                        fieldnames_CONST  = fieldnames(new_class.CONST); %list of field in CONST
                        for ii = 1:size(fieldnames_CONST,1)
                            if isfield(provider.CONST, fieldnames_CONST{ii,1})
                                new_class.CONST.(fieldnames_CONST{ii,1}) = provider.CONST.(fieldnames_CONST{ii,1});
                            %else
                                %warning(['Constant "' fieldnames_CONST{ii,1} '" in class "' class(new_class) '" not populated.'])
                            end
                        end
                    end
                    
                    %provide STATVAR, will be initialized later
                    new_class = provide_STATVAR(new_class);
                    
                    %initialize PARA
                    new_class = provide_PARA(new_class);
                    if ~isempty(new_class.PARA)
                        fieldnames_PARA  = fieldnames(new_class.PARA); %list of field in PARA, find the field and populate it
                        for ii = 1:size(fieldnames_PARA,1)
                            fieldname = fieldnames_PARA{ii};
                            if ~isfield(section{k}, fieldname)
                                %warning(['Parameter "' fieldname '" in class "' class(new_class) '" not populated.'])
                                continue
                            end
                            
                            contents = section{k}.(fieldname);
                            
                            if iscell(contents)   
                                % it is a cell array
                                % which means it is some sort of list or 1D array
                                if ~isempty(contents) && isnumeric(contents{1})
                                    % it is a number array
                                    new_class.PARA.(fieldname) = cell2mat(contents)';
                                else
                                    % it is a cell array of strings
                                    new_class.PARA.(fieldname) = contents';
                                end
                            elseif isstruct(contents)
                                % it is a structure
                                % which means it is some sort of composit
                                % data type, like a STRAT_MATRIX, V_MATRIX
                                % or similar.
                                if isfield(contents,'type') 
                                    if strcmp(contents.type,'STRAT_MATRIX')
                                        % it is a STRAT_MATRIX, a table
                                        % with columnar data, which must
                                        % have 'depth' in the first column
                                        
                                        if ~strcmp(contents.names{1},'depth')
                                            error(['STRAT_MATRIX does not have "depth" as first column. "' fieldname '" in class "' class(new_class) '" not populated.'])
                                        end
                                        
                                        % If matrix contains any string values,
                                        % it is read as a 1D cell array of 1D cell arrays.
                                        if size(contents.values,1) == 1 && iscell(contents.values{1})                                
                                            % this matrix contains some text data
                                            % handle this special case
                                            contents.values = vertcat(contents.values{:});
                                            % now it has same format as numeric matrix
                                        end
                                        
                                        for kk = 1:length(contents.names)
                                            cname = contents.names{kk};

                                            % check type, just to be sure
                                            if isnumeric(contents.values{1,kk})
                                                % it is a numeric type, so convert
                                                % to matrix
                                                new_class.PARA.(fieldname).(cname) = cell2mat(contents.values(:,kk));
                                            else
                                                % it contains text, keep as cell array
                                                new_class.PARA.(fieldname).(cname) = contents.values(:,kk);
                                            end
                                        end
                                    elseif strcmp(contents.type,'V_MATRIX')
                                        % it is a V_MATRIX, a table with columnar data

                                        % If matrix contains any string values,
                                        % it is read as a 1D cell array of 1D cell arrays.
                                        if size(contents.values,1) == 1 && iscell(contents.values{1})                                
                                            % this matrix contains some text data
                                            % handle this special case
                                            contents.values = vertcat(contents.values{:});
                                            % now it has same format as numeric matrix
                                        end

                                        for kk = 1:length(contents.names)
                                            cname = contents.names{kk};

                                            % check type, just to be sure
                                            if isnumeric(contents.values{1,kk})
                                                % it is a numeric type, so convert
                                                % to matrix
                                                new_class.PARA.(fieldname).(cname) = cell2mat(contents.values(:,kk));
                                            else
                                                % it contains text, keep as cell array
                                                new_class.PARA.(fieldname).(cname) = contents.values(:,kk);
                                            end
                                        end
                                    elseif strcmp(contents.type,'MATRIX')
                                        % it is a MATRIX, a table with columnar data
                                        % but it has no headers/column
                                        % names

                                        % NB THIS CODE IS NOT TESTED!!!!

                                        % If matrix contains any string values,
                                        % it is read as a 1D cell array of 1D cell arrays.
                                        if size(contents.values,1) == 1 && iscell(contents.values{1})                                
                                            % this matrix contains some text data
                                            % handle this special case
                                            contents.values = vertcat(contents.values{:});
                                            % now it has same format as numeric matrix
                                        end
                                        
                                        if isnumeric(contents.values{1,1})
                                            new_class.PARA.(fieldname) = cell2mat(contents.values);
                                        else
                                            new_class.PARA.(fieldname) = contents.values;
                                        end
                                    else
                                        warning(['Parameter type "' contents.type '" not recognized. Parameter "' fieldname '" in class "' class(new_class) '" not populated.'])
                                    end
                                else
                                    warning(['Unrecognized compound parameter format. Parameter "' fieldname '" in class "' class(new_class) '" not populated.'])
                                end
                            elseif ismatrix(contents) && isempty(contents)
                                new_class.PARA.(fieldname) = [];
                            elseif isempty(contents)
                                new_class.PARA.(fieldname) = NaN;
                            else
                                % It is just a value or string, keep it
                                new_class.PARA.(fieldname) = contents;
                            end
                        end
                    end
                    provider.CLASSES.(class_name){class_index,1} = new_class;
                    
                    if strcmp(class_sections{i}, 'RUN_INFO') && class_index == 1
                        provider.RUN_INFO_CLASS = new_class;
                    end
                end
            end
        end
        
        %parallel runs
        function provider = update_parameter_file_yaml(provider, param_file_number)
            %change parameter file
            parameter_file = [provider.PARA.result_path provider.PARA.run_name '/' provider.PARA.run_name '_' num2str(param_file_number) '.yml'];
            provider.PARA.parameter_file = parameter_file;
        end
        
        
        function provider = update_run_name_yaml(provider, worker_number)
            %change run name and make new result directories if necessary
            provider.PARA.run_name = [provider.PARA.run_name '_' num2str(worker_number)];
            if ~(exist([provider.PARA.result_path provider.PARA.run_name ], 'dir')==7)
                mkdir([provider.PARA.result_path provider.PARA.run_name]);
            end
        end
        
    end
end

