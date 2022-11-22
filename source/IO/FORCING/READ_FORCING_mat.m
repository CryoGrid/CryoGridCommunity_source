classdef READ_FORCING_mat < READ_FORCING_base
    
    properties
        
    end
    
    methods
        function [data, times] = read_mat(full_file, variables)    %, start_date, end_date)
            
            temp = load(full_file, 'FORCING');
   
            % Implement some cropping of data if start_date and end_date
            % are passed

            for i=1:size(variables,1)
                if isfield(temp, variables{i,1})
                    data.(variables{i,1}) = squeeze(temp.(variables{i,1}));
                end
            end

            times = data.t_span;
        end

    end    
    
end
