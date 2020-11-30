classdef PPROVIDER_EXCEL
    
    properties
        PARA
        CONST
        CLASSES %SUBSURFACE_CLASS, LATERAL_IA_CLASS for which there can be more than one per tile     
        RUN_INFO_CLASS
    end
    
    methods
        
        function pprovider = PPROVIDER_EXCEL(run_name, result_path, constant_file, forcing_path)
            %read constants in CONST filed = all available constants
%             pprovider.PARA.constant_file_full_path = [result_path run_name '/' constant_file '.xlsx'];
%             pprovider.PARA.parameter_file_full_path = [result_path run_name '/' run_name '.xlsx'];
            
            constant_file = [result_path run_name '/' constant_file '.xlsx'];
            parameter_file = [result_path run_name '/' run_name '.xlsx'];

            pprovider.PARA.run_name = run_name;
            pprovider.PARA.result_path = result_path;
            pprovider.PARA.parameter_file = parameter_file;
            pprovider.PARA.constant_file = constant_file;
            pprovider.PARA.forcing_path = forcing_path;
        end
       
        
        function pprovider = read_const(pprovider)
            
            data = read_excel2cell(pprovider, pprovider.PARA.constant_file);
            for i=1:size(data,1)
                if ~isnan(data{i,1}) & ~isnan(data{i,2})
                    pprovider.CONST.(data{i,1}) = data{i,2};
                end
            end
        end
        
        function pprovider = read_parameters(pprovider)
            
            %read the parameter file
            data = read_excel2cell(pprovider, pprovider.PARA.parameter_file);
            i=1;
            while i<=size(data,1)
                if strcmp(data{i,2}, 'index') %class detected
                    functional_group = data{i,1};
                    class_name = data{i+1,1};
                    class_index = data{i+1,2};
                    class_func = str2func(class_name);
                    new_class = class_func();  %initialize the class %CHANGE THIS TO DEFAULT CONSTRUCTOR
                    
                    %initialize CONST
                    new_class = provide_CONST(new_class);
                    if ~isempty(new_class.CONST)
                        fieldnames_CONST  = fieldnames(new_class.CONST); %list of field in PARA, run through the Excel struct, find the field and populate it
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
                        fieldnames_PARA  = fieldnames(new_class.PARA); %list of field in PARA, run through the Excel struct, find the field and populate it
                        for ii = 1:size(fieldnames_PARA,1)
                            j=i+2;
                            while j<=size(data,1) && ~strcmp(data{j,1}, 'CLASS_END')
                                if strcmp(data{j,1}, fieldnames_PARA{ii,1})
                                    %add new keywords here
                                    if strcmp(data{j,2}, 'H_LIST')  %horizontal list
                                        k = 3;
                                        while ~strcmp(data{j,k}, 'END')
                                            k=k+1;
                                        end
                                        if isnumeric(data{j,3})
                                            new_class.PARA.(fieldnames_PARA{ii,1}) = cell2mat(data(j,3:k-1)');
                                        else
                                            new_class.PARA.(fieldnames_PARA{ii,1}) = data(j,3:k-1)';
                                        end
                                        
                                    elseif strcmp(data{j,2}, 'V_LIST') %vertical list, must be in second column
                                        k = j+1;
                                        while ~strcmp(data{k,2}, 'END')
                                            k=k+1;
                                        end
                                        if isnumeric(data{j+1,2})
                                            new_class.PARA.(fieldnames_PARA{ii,1}) = cell2mat(data(j+1:k-1,2));
                                        else
                                            new_class.PARA.(fieldnames_PARA{ii,1}) = data(j+1:k-1,2);
                                        end
                                    elseif strcmp(data{j,2}, 'STRAT_MATRIX') %stratigraphy matrix, read as column vectors, depth in second column of Excel file
                                        k = j+1;
                                        while ~strcmp(data{k,2}, 'END')
                                            k=k+1;
                                        end
                                        l = 3;
                                        while ~strcmp(data{j,l}, 'END')
                                            l=l+1;
                                        end
                                        new_class.PARA.(fieldnames_PARA{ii,1}).depth = cell2mat(data(j+1:k-1,2));
                                        for jj=3:l-1
                                            if isnumeric(data{j+1,jj})
                                                new_class.PARA.(fieldnames_PARA{ii,1}).(data{j,jj}) = cell2mat(data(j+1:k-1,jj));
                                            else
                                                new_class.PARA.(fieldnames_PARA{ii,1}).(data{j,jj}) = data(j+1:k-1,jj);
                                            end
                                        end
                                    elseif strcmp(data{j,2}, 'V_MATRIX') %matrix read as column vectors
                                        k = j+1;
                                        while ~strcmp(data{k,2}, 'END')
                                            k=k+1;
                                        end
                                        l = 3;
                                        while ~strcmp(data{j,l}, 'END')
                                            l=l+1;
                                        end
                                        for jj=3:l-1
                                            if isnumeric(data{j+1,jj})
                                                new_class.PARA.(fieldnames_PARA{ii,1}).(data{j,jj}) = cell2mat(data(j+1:k-1,jj));
                                            else
                                                new_class.PARA.(fieldnames_PARA{ii,1}).(data{j,jj}) = data(j+1:k-1,jj);
                                            end
                                        end
                                    elseif strcmp(data{j,2}, 'MATRIX') %matrix read as column vectors
                                        k = j+1;
                                        while ~strcmp(data{k,2}, 'END')
                                            k=k+1;
                                        end
                                        l = 3;
                                        while ~strcmp(data{j,l}, 'END')
                                            l=l+1;
                                        end

                                        if isnumeric(data{j+1,3})
                                            new_class.PARA.(fieldnames_PARA{ii,1}) = cell2mat(data(j+1:k-1,3:l-1));
                                        else
                                            new_class.PARA.(fieldnames_PARA{ii,1}) = data(j+1:k-1,3:l-1);
                                        end

                                    else  %"normal" field
                                        new_class.PARA.(fieldnames_PARA{ii,1}) = data{j,2};
                                    end
                                    break
                                end
                                j=j+1;
                            end
                            if strcmp(data{j,1}, 'CLASS_END')
                                disp(['WARNING: ' fieldnames_PARA{ii,1} ' not initialized'])
                            end
                        end
                        
                    end
                    %mandatory function in all classes compatible with
                    %PPROVIDER_EXCEL
                    new_class = initialize_excel(new_class);
                    
%                     if strcmp(data{i,2}, 'class_index')
                    pprovider.CLASSES.(class_name){class_index,1} = new_class;
                    if strcmp(functional_group, 'RUN_INFO') && class_index == 1
                        pprovider.RUN_INFO_CLASS = new_class;
                    end

%                     elseif strcmp(data{i,2}, 'function_index')
%                         pprovider.CLASSES.(functional_group){class_index,1} = new_class;
%                     end
                    
                    while i<=size(data,1) && ~strcmp(data{i,1}, 'CLASS_END')
                        i=i+1;
                    end
                end
                
                i=i+1;
            end
            
        end
        
        function [run_info, pprovider] = run(pprovider)
            run_info = copy(pprovider.RUN_INFO_CLASS);
            
            run_info.PPROVIDER = pprovider;
        end
        
        %parallel runs
        function pprovider = update_parameter_file(pprovider, param_file_number)
            %change parameter file
            parameter_file = [pprovider.PARA.result_path pprovider.PARA.run_name '/' pprovider.PARA.run_name '_' num2str(param_file_number) '.xlsx'];
            pprovider.PARA.parameter_file = parameter_file;
        end
        
        function pprovider = update_run_name(pprovider, worker_number)
            %change run name and make new result directories if necessary
            pprovider.PARA.run_name = [pprovider.PARA.run_name '_' num2str(worker_number)];
            if ~(exist([pprovider.PARA.result_path pprovider.PARA.run_name ])==7)
                mkdir([pprovider.PARA.result_path pprovider.PARA.run_name]);
            end
        end
        
        %service functions
        function result = read_excel2cell(pprovider, filename)
            [dummy1, dummy2, result] = xlsread(filename);
        end
        
        
    end
end

