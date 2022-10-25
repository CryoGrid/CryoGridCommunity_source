classdef PROVIDER_EXCEL_edu < PROVIDER_EXCEL
    
    
    methods
        
        function provider = assign_paths_excel_edu(provider, run_name, result_path, constant_file)

            constant_file = [result_path run_name '/' constant_file '.xlsx'];
            parameter_file = [result_path run_name '/' run_name '.xlsx'];
            parameter_file_edu = [result_path run_name '/' run_name '_edu.xlsx'];

            provider.PARA.run_name = run_name;
            provider.PARA.result_path = result_path;
            provider.PARA.parameter_file = parameter_file;
            provider.PARA.parameter_file_edu = parameter_file_edu;
            provider.PARA.constant_file = constant_file;
            %provider.PARA.forcing_path = forcing_path;
        end

        
        function provider = read_const_excel_edu(provider)
            provider = read_const_excel(provider);
        end
        
        function provider = read_parameters_excel_edu(provider)
            provider = read_parameters_excel(provider);
            %read adn overwrite
            data = read_excel2cell(provider, provider.PARA.parameter_file_edu);
            disp(provider.PARA.parameter_file_edu)
            for i=1:size(data,1)
                class_name = data{i,1};
                class_index = data{i,2};
                PARA_name = data{i,3};
                PARA_value = data{i,4};
                
                try
%                     if isnumeric(PARA_value)
%                         provider.CLASSES.(class_name){class_index,1}.(PARA_name) = PARA_value;
%                         new_class.PARA.(fieldnames_PARA{ii,1}) = cell2mat(data(j+1:k-1,2));
%                     else
%                         new_class.PARA.(fieldnames_PARA{ii,1}) = data(j+1:k-1,2);
%                     end
                    provider.CLASSES.(class_name){class_index,1}.PARA.(PARA_name) = PARA_value;
                end
            end
        end
        
        
%         function provider = read_parameters_excel(provider)
%             
%             %read the parameter file
%             data = read_excel2cell(provider, provider.PARA.parameter_file);
%             disp(provider.PARA.parameter_file)
%             i=1;
%             while i<=size(data,1)
%                 if strcmp(data{i,2}, 'index') %class detected
%                     functional_group = data{i,1};
%                     class_name = data{i+1,1};
%                     class_index = data{i+1,2};
%                     class_func = str2func(class_name);
%                     new_class = class_func();  %initialize the class %CHANGE THIS TO DEFAULT CONSTRUCTOR
%                     
%                     %initialize CONST
%                     new_class = provide_CONST(new_class);
%                     if ~isempty(new_class.CONST)
%                         fieldnames_CONST  = fieldnames(new_class.CONST); %list of field in PARA, run through the Excel struct, find the field and populate it
%                         for ii = 1:size(fieldnames_CONST,1)
%                             if isfield(provider.CONST, fieldnames_CONST{ii,1})
%                                 new_class.CONST.(fieldnames_CONST{ii,1}) = provider.CONST.(fieldnames_CONST{ii,1});
%                             else
%                                % disp(['WARNING: constant ' fieldnames_CONST{ii,1} ' in class ' class(new_class) ' not populated.'])
%                             end
%                         end
%                     end
%                     
%                     %provide STATVAR, will be initialized later
%                     new_class = provide_STATVAR(new_class);
%                     
%                     %initialize PARA
%                     new_class = provide_PARA(new_class);
%                     if ~isempty(new_class.PARA)
%                         fieldnames_PARA  = fieldnames(new_class.PARA); %list of fields in the PARA structure of new_class
%                         
%                         % Iterate over all fieldnames in new_class.PARA
%                         for ii = 1:size(fieldnames_PARA,1)
%                             %run through the Excel struct, find the field and populate it
%                             j=i+2;
%                             while j<=size(data,1) && ~strcmp(data{j,1}, 'CLASS_END')
%                                 if strcmp(data{j,1}, fieldnames_PARA{ii,1})
%                                     %add new keywords here
%                                     if strcmp(data{j,2}, 'H_LIST')  %horizontal list
%                                         k = 3;
%                                         while ~strcmp(data{j,k}, 'END')
%                                             k=k+1;
%                                         end
%                                         if isnumeric(data{j,3})
%                                             new_class.PARA.(fieldnames_PARA{ii,1}) = cell2mat(data(j,3:k-1)');
%                                         else
%                                             new_class.PARA.(fieldnames_PARA{ii,1}) = data(j,3:k-1)';
%                                         end
%                                         
%                                     elseif strcmp(data{j,2}, 'V_LIST') %vertical list, must be in second column
%                                         k = j+1;
%                                         while ~strcmp(data{k,2}, 'END')
%                                             k=k+1;
%                                         end
%                                         if isnumeric(data{j+1,2})
%                                             new_class.PARA.(fieldnames_PARA{ii,1}) = cell2mat(data(j+1:k-1,2));
%                                         else
%                                             new_class.PARA.(fieldnames_PARA{ii,1}) = data(j+1:k-1,2);
%                                         end
%                                     elseif strcmp(data{j,2}, 'STRAT_MATRIX') %stratigraphy matrix, read as column vectors, depth in second column of Excel file
%                                         k = j+1;
%                                         while ~strcmp(data{k,2}, 'END')
%                                             k=k+1;
%                                         end
%                                         l = 3;
%                                         while ~strcmp(data{j,l}, 'END')
%                                             l=l+1;
%                                         end
%                                         new_class.PARA.(fieldnames_PARA{ii,1}).depth = cell2mat(data(j+1:k-1,2));
%                                         for jj=3:l-1
%                                             if isnumeric(data{j+1,jj})
%                                                 new_class.PARA.(fieldnames_PARA{ii,1}).(data{j,jj}) = cell2mat(data(j+1:k-1,jj));
%                                             else
%                                                 new_class.PARA.(fieldnames_PARA{ii,1}).(data{j,jj}) = data(j+1:k-1,jj);
%                                             end
%                                         end
%                                     elseif strcmp(data{j,2}, 'V_MATRIX') %matrix read as column vectors
%                                         k = j+1;
%                                         while ~strcmp(data{k,2}, 'END')
%                                             k=k+1;
%                                         end
%                                         l = 3;
%                                         while ~strcmp(data{j,l}, 'END')
%                                             l=l+1;
%                                         end
%                                         for jj=3:l-1
%                                             if isnumeric(data{j+1,jj})
%                                                 new_class.PARA.(fieldnames_PARA{ii,1}).(data{j,jj}) = cell2mat(data(j+1:k-1,jj));
%                                             else
%                                                 new_class.PARA.(fieldnames_PARA{ii,1}).(data{j,jj}) = data(j+1:k-1,jj);
%                                             end
%                                         end
%                                     elseif strcmp(data{j,2}, 'MATRIX') %matrix read as column vectors
%                                         k = j+1;
%                                         while ~strcmp(data{k,2}, 'END')
%                                             k=k+1;
%                                         end
%                                         l = 3;
%                                         while ~strcmp(data{j,l}, 'END')
%                                             l=l+1;
%                                         end
% 
%                                         if isnumeric(data{j+1,3})
%                                             new_class.PARA.(fieldnames_PARA{ii,1}) = cell2mat(data(j+1:k-1,3:l-1));
%                                         else
%                                             new_class.PARA.(fieldnames_PARA{ii,1}) = data(j+1:k-1,3:l-1);
%                                         end
% 
%                                     else  %"normal" field
%                                         new_class.PARA.(fieldnames_PARA{ii,1}) = data{j,2};
%                                     end
%                                     break
%                                 end
%                                 j=j+1;
%                             end
% %                             if strcmp(data{j,1}, 'CLASS_END')
% %                                 disp(['WARNING: ' fieldnames_PARA{ii,1} ' not initialized'])
% %                             end
%                         end
%                         
%                     end
%                     %mandatory function in all classes compatible with
%                     %PROVIDER_EXCEL - REMOVED!!! - no longer necessary
%                     %new_class = initialize_excel(new_class);
%                     
%                     provider.CLASSES.(class_name){class_index,1} = new_class;
%                     if strcmp(functional_group, 'RUN_INFO') && class_index == 1
%                         provider.RUN_INFO_CLASS = new_class;
%                     end
%                     
%                     while i<=size(data,1) && ~strcmp(data{i,1}, 'CLASS_END')
%                         i=i+1;
%                     end
%                 end
%                 
%                 i=i+1;
%             end
%             
%         end
        
        
        %parallel runs
%         function provider = update_parameter_file_excel_edu(provider, param_file_number)
%             %change parameter file
%             parameter_file = [provider.PARA.result_path provider.PARA.run_name '/' provider.PARA.run_name '_' num2str(param_file_number) '.xlsx'];
%             provider.PARA.parameter_file = parameter_file;
%         end
%         
%         function provider = update_run_name_excel_edu(provider, worker_number)
%             %change run name and make new result directories if necessary
%             provider.PARA.run_name = [provider.PARA.run_name '_' num2str(worker_number)];
%             if ~(exist([provider.PARA.result_path provider.PARA.run_name ])==7)
%                 mkdir([provider.PARA.result_path provider.PARA.run_name]);
%             end
%         end
        
        %service functions
%         function result = read_excel2cell(provider, filename)
%             [dummy1, dummy2, result] = xlsread(filename);
%         end
        
        
    end
end

