classdef READ_FORCING_mat < READ_FORCING_base
    
    properties
        
    end
    
    methods
        function [data, timestamp] = read_mat(forcing, full_file, variables)    %, start_date, end_date)
            
            temp = load(full_file, 'FORCING');
            temp=temp.FORCING.data;
            % Implement some cropping of data if start_date and end_date
            % are passed

            for i=1:size(variables,1)
                if isfield(temp, variables{i,1})
                    data.(variables{i,1}) = squeeze(temp.(variables{i,1}));
                end
            end
            if isfield(temp, 't_span')
                timestamp = temp.t_span;
            elseif isfield(temp, 'timeForcing')
               timestamp = temp.timeForcing;
            else
                timestamp = 0;
                disp('no suitable timestamp');
            end
        end
        
%         %for TopoScale
        function data = read_mat_ERA4D(forcing, full_file)
            data = load(full_file, 'era');
            data = data.era;
        end


    end    
    
end
