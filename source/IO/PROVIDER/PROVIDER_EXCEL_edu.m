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

        
    end
end

