classdef PPROVIDER_YAML
    
    properties
        PARA
        CONST
        CLASSES %SUBSURFACE_CLASS, LATERAL_IA_CLASS for which there can be more than one per tile     
        RUN_INFO_CLASS
    end
    
    methods
        
        function pprovider = PPROVIDER_YAML(run_name, result_path, constant_file, forcing_path)
            
            constant_file = [result_path run_name '/' constant_file '.yml'];
            parameter_file = [result_path run_name '/' run_name '.yml'];

            pprovider.PARA.run_name = run_name;
            pprovider.PARA.result_path = result_path;
            pprovider.PARA.parameter_file = parameter_file;
            pprovider.PARA.constant_file = constant_file;
            pprovider.PARA.forcing_path = forcing_path;
        end
       
        
        function pprovider = read_const(pprovider)
            
            data = yml.ReadYaml(pprovider.PARA.constant_file);
            fnames = fieldnames(data.CONST);
            for i=1:length(fieldnames)
                if ~isnan(data.CONST.(fnames{i}))
                    pprovider.CONST.(fnames{i}) = data.CONST.(fnames{i});
                end
            end
        end
        
        function pprovider = read_parameters(pprovider)
            
            %read the parameter file
            data = yml.ReadYaml(pprovider.PARA.parameter_file);

            % get all the CLASS sections defined in the file
            class_sections = fieldnames(data); 
            
            % loop over all sections
            for i = 1:length(class_sections)
                section = class_sections{1};
                
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
                            if isfield(pprovider.CONST, fieldnames_CONST{ii,1})
                                new_class.CONST.(fieldnames_CONST{ii,1}) = pprovider.CONST.(fieldnames_CONST{ii,1});
                            else
                                disp(['WARNING: constant ' fieldnames_CONST{ii,1} ' in class ' class(new_class) ' not populated.'])
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
                            if ~isfield(pprovider.PARA, fieldname)
                                disp(['WARNING: parameter ' fieldname ' in class ' class(new_class) ' not populated.'])
                                continue
                            end
                            
                            contents = section{k}.(fieldname);
                            
                            if iscell(contents)   
                                % it is a cell array
                                % which means it is some sort of list or 1D array
                                if ~isempty(contents) && isnumeric(contents{1})
                                    % it is a number array
                                    new_class.PARA.(fieldname) = cell2mat(contents);
                                else
                                    % it is a cell array of strings
                                    new_class.PARA.(fieldname) = contents;
                                end
                            elseif isstruct(contents)
                                % it is a structure
                                % which means it is a table with columnar data
                                for kk = 1:length(contents.names)
                                    cname = contents.names{kk};
                                    if isnumeric(contents.values{1,kk})
                                        % it is a numeric type, so convert
                                        % to matrix
                                        new_class.PARA.(fieldname).(cname) = cell2mat(contents.values(:,kk));
                                    else
                                        % it contains text, keep as cell
                                        % array
                                        new_class.PARA.(fieldname).(cname) = contents.values(:,kk);
                                    end
                                end
                            else
                                % It is just a value or string, keep it
                                new_class.PARA.(fieldname) = contents;
                            end
                        end
                    end
                    %mandatory function in all classes compatible with
                    %PPROVIDER_YAML
                    new_class = initialize_yaml(new_class);

                    pprovider.CLASSES.(class_name){class_index,1} = new_class;
                    
                    if strcmp(section, 'RUN_INFO') && class_index == 1
                        pprovider.RUN_INFO_CLASS = new_class;
                    end
                end
            end
        end

        
        function [run_info, pprovider] = run(pprovider)
            run_info = copy(pprovider.RUN_INFO_CLASS);
            run_info.PPROVIDER = pprovider;
        end
        
        
        %parallel runs
        function pprovider = update_parameter_file(pprovider, param_file_number)
            %change parameter file
            parameter_file = [pprovider.PARA.result_path pprovider.PARA.run_name '/' pprovider.PARA.run_name '_' num2str(param_file_number) '.yml'];
            pprovider.PARA.parameter_file = parameter_file;
        end
        
        
        function pprovider = update_run_name(pprovider, worker_number)
            %change run name and make new result directories if necessary
            pprovider.PARA.run_name = [pprovider.PARA.run_name '_' num2str(worker_number)];
            if ~(exist([pprovider.PARA.result_path pprovider.PARA.run_name ])==7)
                mkdir([pprovider.PARA.result_path pprovider.PARA.run_name]);
            end
        end
        
    end
end

