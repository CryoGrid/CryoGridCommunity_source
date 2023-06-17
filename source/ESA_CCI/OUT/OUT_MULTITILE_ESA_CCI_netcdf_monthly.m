classdef OUT_MULTITILE_ESA_CCI_netcdf_monthly < matlab.mixin.Copyable
    %only works for yearly output

    properties
        
        PARA
        STATVAR
        CONST
        TEMP
        OUTPUT_TIME
        SAVE_TIME
    end
    
    methods
        
        
        function out = provide_PARA(out)
            out.PARA.out_folder = [];
            out.PARA.output_timestep_months = []; %in months
        end
        
        function out = provide_CONST(out)
            %out.CONST.day_sec = [];
        end
        
        function out = provide_STATVAR(out)

        end
        
        
        function out = finalize_init(out, tile)

            if ~(exist([out.PARA.out_folder tile.RUN_INFO.PARA.run_name])==7)
                mkdir([out.PARA.out_folder tile.RUN_INFO.PARA.run_name]);
            end
            out.PARA.out_folder = [out.PARA.out_folder tile.RUN_INFO.PARA.run_name '/'];
            
            out = reset_STATVAR(out); %initializes STATVAR
            
            variable_names = fieldnames(out.STATVAR);
            model_years = [year(tile.FORCING.PARA.start_time):year(tile.FORCING.PARA.end_time)]';
            
            target_dates=[];
            for month=1:out.PARA.output_timestep_months:12
                target_dates = [target_dates datenum(model_years, month, 1)];
            end
            target_dates=target_dates';
            target_dates=target_dates(:);
            
            target_dates(target_dates < tile.FORCING.PARA.start_time | target_dates > tile.FORCING.PARA.end_time) = [];
            
            if exist([out.PARA.out_folder 'out_monthly_' num2str(tile.PARA.range(1,1)) '_' num2str(tile.PARA.range(end,1)) '.nc'])==2
                delete([out.PARA.out_folder 'out_monthly_' num2str(tile.PARA.range(1,1)) '_' num2str(tile.PARA.range(end,1)) '.nc'])
            end
            for i=1:size(variable_names,1)
                nccreate([out.PARA.out_folder 'out_monthly_' num2str(tile.PARA.range(1,1)) '_' num2str(tile.PARA.range(end,1)) '.nc'], variable_names{i,1}, 'Format', 'netcdf4', 'Datatype', 'uint16', 'Dimensions', ...
                    {'x',size(tile.PARA.range,1) ,'y', size(target_dates,1), 'z', tile.PARA.ensemble_size}, 'FillValue', 2^16-1);
            end
            
            out.SAVE_TIME = [target_dates; Inf];
            out.TEMP.date_count = 1;
        end
        
        function out = store_OUT(out, tile) 

            if tile.t>=out.SAVE_TIME(1,1)

                out.STATVAR.current_T1m = tile.SUBSURFACE_CLASS.STATVAR.T(14,:);

                T_pos = double(tile.SUBSURFACE_CLASS.STATVAR.T(5:end,:)>0);
                T_pos=[zeros(1, size(T_pos,2)) ; T_pos];
                change_index = T_pos(1:end-1,:)-T_pos(2:end,:); %1 for unfrozen->frozen, -1 for frozen->unfrozen %N
                AL_cell = double(cumsum(double(change_index==1),1)==0); %N
                layerThick = tile.SUBSURFACE_CLASS.STATVAR.layerThick(5:end,:);
                layerThick(1,:) = tile.SUBSURFACE_CLASS.STATVAR.layerThick_first_ground_cell;
                out.STATVAR.current_thawDepth = sum(layerThick .* AL_cell, 1) .* double(tile.SUBSURFACE_CLASS.STATVAR.T(5,:)>0);

                variable_names = fieldnames(out.STATVAR);
                for i=1:size(variable_names,1)
                    datapackage = out.STATVAR.(variable_names{i,1});

                    datapackage = reshape(datapackage, size(datapackage, 2) ./ tile.PARA.ensemble_size, 1, tile.PARA.ensemble_size);

                    datapackage = uint16( round( (datapackage - out.TEMP.scale_offset(i,1)) ./ out.TEMP.scale_factor(i,1) ) );

                    ncwrite([out.PARA.out_folder 'out_monthly_' num2str(tile.PARA.range(1,1)) '_' num2str(tile.PARA.range(end,1)) '.nc'], ...
                        variable_names{i,1}, datapackage, [1 out.TEMP.date_count 1], [1 1 1]);

                end

                out.TEMP.date_count = out.TEMP.date_count + 1;


                out = reset_STATVAR(out);
                out.SAVE_TIME(1,:) = [];

            end


        end
        

            
        function out = reset_STATVAR(out)
            out.STATVAR.current_T1m = [];
            out.STATVAR.current_thawDepth  = [];

            out.TEMP.scale_factor = [0.002; 0.01];
            out.TEMP.scale_offset = [-70; 0];
            
          %  out.TEMP.scale_factor = [0.002;0.002;0.002; 0.002;0.002;0.002; 0.01; 1; 1; 0.01; 1; 0.002;0.002;0.002;0.002;0.002;0.002];
          %  out.TEMP.scale_offset = [-70;-70;-70;-70;-70;-70; 0; -10; 0; 0; 0; -70;-70;-70;-70;-70;-70];
        end

    end
end

