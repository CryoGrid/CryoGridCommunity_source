classdef PROVIDER_MAT < BASE_PROVIDER
    
    
    methods
        
        function provider = assign_paths_mat(provider, run_name, result_path)

            
            parameter_file = [result_path run_name '/' run_name '.mat'];
            provider.PARA.run_name = run_name;
            provider.PARA.result_path = result_path;
            provider.PARA.parameter_file = parameter_file;
            
            provider_from_file = load(provider.PARA.parameter_file);

            fn = fieldnames(provider_from_file);
            provider_from_file = provider_from_file.(fn{1,1});
            
            fn = fieldnames(provider_from_file);
            
            for i=1:size(fn,1)
               if any(strcmp(properties(provider), fn{i,1})) 
                  fn2 = fieldnames(provider_from_file.(fn{i,1}));
                  for j=1:size(fn2,1)
                      provider.(fn{i,1}).(fn2{j,1}) = provider_from_file.(fn{i,1}).(fn2{j,1});
                  end
               end
            end

        end
    end
end

