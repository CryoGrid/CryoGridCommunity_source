classdef OUT_MULTITILE_ESA_CCI_netcdf < matlab.mixin.Copyable
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
            out.PARA.output_timestep = []; %in days, if empty, average over entire period
            out.PARA.save_date = [];
            out.PARA.save_interval = []; %in years, this is when the output file is written
        end
        
        function out = provide_CONST(out)
            out.CONST.day_sec = [];
        end
        
        function out = provide_STATVAR(out)

        end
        
        
        function out = finalize_init(out, tile)
% 
%                       tile.RUN_INFO.PARA.run_name = '';
%                       tile.PARA.range =[1; 2];

            if ~(exist([out.PARA.out_folder tile.RUN_INFO.PARA.run_name])==7)
                mkdir([out.PARA.out_folder tile.RUN_INFO.PARA.run_name]);
            end
            out.PARA.out_folder = [out.PARA.out_folder tile.RUN_INFO.PARA.run_name '/'];
            
            out = set2zero(out);
            out = reset_STATVAR(out); %initializes STATVAR
            
            variable_names = fieldnames(out.STATVAR);
            number_of_years = str2num(datestr(tile.FORCING.PARA.end_time, 'yyyy')) - str2num(datestr(tile.FORCING.PARA.start_time, 'yyyy'))+1;
            
            if exist([out.PARA.out_folder 'out_' num2str(tile.PARA.range(1,1)) '_' num2str(tile.PARA.range(end,1)) '.nc'])==2
                delete([out.PARA.out_folder 'out_' num2str(tile.PARA.range(1,1)) '_' num2str(tile.PARA.range(end,1)) '.nc'])
            end
            for i=1:size(variable_names,1)
                nccreate([out.PARA.out_folder 'out_' num2str(tile.PARA.range(1,1)) '_' num2str(tile.PARA.range(end,1)) '.nc'], variable_names{i,1}, 'Format', 'netcdf4', 'Datatype', 'uint16', 'Dimensions', ...
                    {'x',size(tile.PARA.range,1) ,'y', number_of_years, 'z', tile.PARA.ensemble_size}, 'FillValue', 2^16-1);
                
            end
            

            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = tile.FORCING.PARA.end_time;
            else
                if datenum([out.PARA.save_date datestr(tile.FORCING.PARA.start_time,'yyyy')], 'dd.mm.yyyy') <= tile.FORCING.PARA.start_time + 30
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(tile.FORCING.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                else
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date datestr(tile.FORCING.PARA.start_time,'yyyy')], 'dd.mm.yyyy'));
                end
            end
            
            if isempty(out.PARA.output_timestep) || isnan(out.PARA.output_timestep)
                out.OUTPUT_TIME = out.SAVE_TIME;
            else
                out.OUTPUT_TIME = tile.FORCING.PARA.start_time + out.PARA.output_timestep;
            end
            
            out.TEMP.year_count = 1;

            tic
        end
        
        function out = store_OUT(out, tile) 

            out.TEMP.acc_GST = out.TEMP.acc_GST + tile.SUBSURFACE_CLASS.STATVAR.T(5,:);
            out.TEMP.acc_T1m = out.TEMP.acc_T1m + tile.SUBSURFACE_CLASS.STATVAR.T(14,:);
            out.TEMP.acc_T2m = out.TEMP.acc_T2m + tile.SUBSURFACE_CLASS.STATVAR.T(22,:);
            out.TEMP.acc_T5m = out.TEMP.acc_T5m + tile.SUBSURFACE_CLASS.STATVAR.T(31,:);
            out.TEMP.acc_T10m = out.TEMP.acc_T10m + tile.SUBSURFACE_CLASS.STATVAR.T(35,:);
            out.TEMP.acc_T20m = out.TEMP.acc_T20m + tile.SUBSURFACE_CLASS.STATVAR.T(37,:);

            %out.TEMP.acc_FT_state = out.TEMP.acc_FT_state + double(tile.SUBSURFACE_CLASS.STATVAR.T(5:35,:) <0); %only check until 10m for talik
            %out.TEMP.acc_FT_isothermal = out.TEMP.acc_FT_isothermal & (tile.SUBSURFACE_CLASS.STATVAR.T(5:35,:) <0 & tile.SUBSURFACE_CLASS.STATVAR.T(5:35,:) >=-0.1); %cells that are always isothermal around 0 degree C
            
            out.TEMP.acc_N = out.TEMP.acc_N + 1;
            
            out.TEMP.acc_T1m_max = max(out.TEMP.acc_T1m_max, tile.SUBSURFACE_CLASS.STATVAR.T(14,:));
            out.TEMP.acc_T1m_min = min(out.TEMP.acc_T1m_min, tile.SUBSURFACE_CLASS.STATVAR.T(14,:));
            out.TEMP.acc_T2m_max = max(out.TEMP.acc_T1m_max, tile.SUBSURFACE_CLASS.STATVAR.T(22,:));
            out.TEMP.acc_T2m_min = min(out.TEMP.acc_T1m_max, tile.SUBSURFACE_CLASS.STATVAR.T(22,:));
            out.TEMP.acc_T1m_square = out.TEMP.acc_T1m_square + tile.SUBSURFACE_CLASS.STATVAR.T(14,:).^2;
            out.TEMP.acc_T2m_square = out.TEMP.acc_T2m_square + tile.SUBSURFACE_CLASS.STATVAR.T(22,:).^2;
            
%             out.TEMP.acc_FDD = out.TEMP.acc_FDD + tile.SUBSURFACE_CLASS.STATVAR.surf_T .* double(tile.SUBSURFACE_CLASS.STATVAR.surf_T<0);  
%             out.TEMP.acc_TDD = out.TEMP.acc_TDD + tile.SUBSURFACE_CLASS.STATVAR.surf_T .* double(tile.SUBSURFACE_CLASS.STATVAR.surf_T>=0);
%            
            T_pos = double(tile.SUBSURFACE_CLASS.STATVAR.T(5:end,:)>0);
            T_pos=[zeros(1, size(T_pos,2)) ; T_pos]; 
            change_index = T_pos(1:end-1,:)-T_pos(2:end,:); %1 for unfrozen->frozen, -1 for frozen->unfrozen %N
            AL_cell = double(cumsum(double(change_index==1),1)==0); %N
%             D = tile.SUBSURFACE_CLASS.STATVAR.K_grid(2:end,:) - profile.K_grid(1:end-1,:); %N
            layerThick = tile.SUBSURFACE_CLASS.STATVAR.layerThick(5:end,:);
            layerThick(1,:) = tile.SUBSURFACE_CLASS.STATVAR.layerThick_first_ground_cell;
            thaw_depth = sum(layerThick .* AL_cell, 1) .* double(tile.SUBSURFACE_CLASS.STATVAR.T(5,:)>0);
            out.TEMP.acc_max_ALT = max(out.TEMP.acc_max_ALT, thaw_depth);
            
            out.TEMP.acc_snow_density = out.TEMP.acc_snow_density + sum(tile.SUBSURFACE_CLASS.STATVAR.ice_snow,1)./max(1e-20, sum(tile.SUBSURFACE_CLASS.STATVAR.layerThick_snow,1)).*917; 
            out.TEMP.acc_snow_depth = out.TEMP.acc_snow_depth + sum(tile.SUBSURFACE_CLASS.STATVAR.layerThick_snow,1);
            out.TEMP.acc_snow_covered_days = out.TEMP.acc_snow_covered_days + double(sum(tile.SUBSURFACE_CLASS.STATVAR.layerThick_snow,1)>0);

             if tile.t >= out.OUTPUT_TIME 
                 
                 out.STATVAR.T1m = [out.STATVAR.T1m; out.TEMP.acc_T1m ./ out.TEMP.acc_N];
                 out.STATVAR.T2m = [out.STATVAR.T2m; out.TEMP.acc_T2m ./ out.TEMP.acc_N];
                 out.STATVAR.T5m = [out.STATVAR.T5m; out.TEMP.acc_T5m ./ out.TEMP.acc_N];
                 out.STATVAR.T10m = [out.STATVAR.T10m; out.TEMP.acc_T10m ./ out.TEMP.acc_N];
                 out.STATVAR.T20m = [out.STATVAR.T20m; out.TEMP.acc_T20m ./ out.TEMP.acc_N];
                 out.STATVAR.GST = [out.STATVAR.GST; out.TEMP.acc_GST ./ out.TEMP.acc_N];
                 out.STATVAR.max_ALT  = [out.STATVAR.max_ALT; out.TEMP.acc_max_ALT];
                 out.STATVAR.snow_density = [out.STATVAR.snow_density; out.TEMP.acc_snow_density ./ out.TEMP.acc_N];
                 out.STATVAR.snow_depth = [out.STATVAR.snow_depth; out.TEMP.acc_snow_depth ./ out.TEMP.acc_N];
                 out.STATVAR.snow_covered_days  = [out.STATVAR.snow_covered_days; out.TEMP.acc_snow_covered_days ./ (out.CONST.day_sec ./ tile.timestep)];
                 out.STATVAR.T1m_max  = [out.STATVAR.T1m_max; out.TEMP.acc_T1m_max];
                 out.STATVAR.T1m_min = [out.STATVAR.T1m_min; out.TEMP.acc_T1m_min];
                 out.STATVAR.T2m_max = [out.STATVAR.T2m_max; out.TEMP.acc_T2m_max];
                 out.STATVAR.T2m_min = [out.STATVAR.T2m_min; out.TEMP.acc_T2m_min];
                 out.STATVAR.T1m_std = [out.STATVAR.T1m_std; real(sqrt(out.TEMP.acc_T1m_square./out.TEMP.acc_N - (out.TEMP.acc_T1m ./ out.TEMP.acc_N).^2))];
                 out.STATVAR.T2m_std = [out.STATVAR.T2m_std; real(sqrt(out.TEMP.acc_T2m_square./out.TEMP.acc_N - (out.TEMP.acc_T2m ./ out.TEMP.acc_N).^2))];
%                  out.STATVAR.FDD = [out.STATVAR.FDD; out.TEMP.acc_FDD]; %initialization
%                  out.STATVAR.TDD = [out.STATVAR.TDD; out.TEMP.acc_FDD];
                 out.STATVAR.FT_state = [out.STATVAR.FT_state; tile.SUBSURFACE_CLASS.STATVAR.FT_state];

                 out = set2zero(out);
                 
                 
                 if tile.t>=out.SAVE_TIME
                                          
                     variable_names = fieldnames(out.STATVAR);
                     for i=1:size(variable_names,1)
                         datapackage = out.STATVAR.(variable_names{i,1});
                         
                         datapackage = reshape(datapackage, size(datapackage, 2) ./ tile.PARA.ensemble_size, 1, tile.PARA.ensemble_size);
                         
                         datapackage = uint16( round( (datapackage - out.TEMP.scale_offset(i,1)) ./ out.TEMP.scale_factor(i,1) ) );

                         ncwrite([out.PARA.out_folder 'out_' num2str(tile.PARA.range(1,1)) '_' num2str(tile.PARA.range(end,1)) '.nc'], ...
                             variable_names{i,1}, datapackage, [1 out.TEMP.year_count 1], [1 1 1]);

                     end
                     
                     out.TEMP.year_count = out.TEMP.year_count + 1;
                     
%                      filename =  [out.PARA.out_folder tile.RUN_INFO.PARA.run_name '/'  'out_' datestr(tile.t-1, 'yyyy') '_' num2str(tile.PARA.range(1,1)) '_' num2str(tile.PARA.range(end,1)) '.mat'];
%                      var = out.STATVAR;
%                      save(filename, 'var');
                     

                     out = reset_STATVAR(out);
                     out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                    
                     out = set2zero(out);
                     toc
                     tic
                 end

                     
                 if isempty(out.PARA.output_timestep) || isnan(out.PARA.output_timestep)
                     out.OUTPUT_TIME = out.SAVE_TIME;
                 else
                     out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                 end

            end

%             out.TEMP.stratigraphy = profile.stratigraphy;
%             out.TEMP.wind_compaction_timescale = profile.wind_compaction_timescale;
%             out.TEMP.sf_factor = profile.snowfall_factor;

        end
        
%         function out = get_yearly_av_OUT(out, year)
%             OUT.year_list = [OUT.year_list; year];
%             OUT.GST_av = [OUT.GST_av; OUT.ACC.GST ./ OUT.ACC.N ];
%             OUT.T1m_av = [OUT.T1m_av; OUT.ACC.T1m ./ OUT.ACC.N ]; 
%             OUT.T2m_av = [OUT.T2m_av; OUT.ACC.T2m ./ OUT.ACC.N ];
%             OUT.T5m_av = [OUT.T5m_av; OUT.ACC.T5m ./ OUT.ACC.N ];
%             OUT.T10m_av = [OUT.T10m_av; OUT.ACC.T10m ./ OUT.ACC.N ];
% 	    OUT.T20m_av = [OUT.T20m_av; OUT.ACC.T20m ./ OUT.ACC.N ];
% 
% 	    OUT.T1m_max = [ OUT.T1m_max; OUT.ACC.T1m_max];
% 	    OUT.T1m_min = [ OUT.T1m_min; OUT.ACC.T1m_min];
% 	    OUT.T2m_max = [ OUT.T2m_max; OUT.ACC.T2m_max];
% 	    OUT.T2m_min = [ OUT.T2m_min; OUT.ACC.T2m_min];
% 
% 	    OUT.T1m_std = [ OUT.T1m_std; sqrt(OUT.ACC.T1m_square./OUT.ACC.N - (OUT.ACC.T1m ./ OUT.ACC.N).^2)];
% 	    OUT.T2m_std = [ OUT.T2m_std; sqrt(OUT.ACC.T2m_square./OUT.ACC.N - (OUT.ACC.T2m ./ OUT.ACC.N).^2)];
% 
%             OUT.ALT_max = [OUT.ALT_max; OUT.ACC.max_ALT];
%             OUT.snow_density =[OUT.snow_density; OUT.ACC.snow_density];
%             OUT.snow_depth = [OUT.snow_depth; OUT.ACC.snow_depth ./ OUT.ACC.snow_covered_days];
%             OUT.snow_covered_days = [OUT.snow_covered_days; OUT.ACC.snow_covered_days];
%             FT_code=OUT.ACC.FT_state ./ OUT.ACC.N;
%             FT_code(FT_code==1)=2;  %frozen
%             FT_code(FT_code>0 & FT_code < 1)=1;  %freeze thaw
% 
% 	    OUT.gain_loose = FT_code.*0;
% 	    OUT.gain_loose(FT_code==2) = 1;  %gain when frozen
% 	    OUT.gain_loose(OUT.ACC.FT_isothermal==1) = 0; %no gain when isothermal
% 	    OUT.gain_loose(FT_code == 1 | FT_code == 0) = -1 ;  %loose when unfrozen or FT, this depends on water table settings for the ensemble member
% 	    for i=1:size(OUT.gain_loose,1)-1  %set cell above AL to gain
% 		OUT.gain_loose(i,:) = OUT.gain_loose(i,:) + double(FT_code(i,:)==1 & FT_code(i+1,:)==2 & ~(OUT.ACC.FT_isothermal(i+1,:)==1)) .*(1 - OUT.gain_loose(i,:)); 
%             end
% 	    
% 
%             FT_code=[ones(1,size(FT_code,2)); FT_code];
%             FT_code = FT_code(1:end-1,:) - FT_code(2:end,:);
% 
%             FT_res = zeros(1, size(FT_code,2));
%             for i = 1: size(FT_code,1)
%                FT_res = FT_res + double (FT_res == 0 & FT_code(i,:) ==-1) .* -1 + double (FT_res == 0 & FT_code(i,:) == 1) .* 1; %initial PF yes no, -1 or 1
%                FT_res = FT_res + double (FT_res >0 & FT_code(i,:) < 0) .* -FT_code(i,:); %switches from no PF, 1, to freeze_thaw or frozen, so 2 or 3 means talik
%             end
%           
%             OUT.FT_state = [OUT.FT_state;  FT_res];
%             
%             dlmwrite([OUT.folder 'GST_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.GST_av])
%             dlmwrite([OUT.folder 'T1m_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T1m_av])
%             dlmwrite([OUT.folder 'T2m_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T2m_av])
%             dlmwrite([OUT.folder 'T5m_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T5m_av])
%             dlmwrite([OUT.folder 'T10m_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T10m_av])
% 	    dlmwrite([OUT.folder 'T20m_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T20m_av])
% 
%             dlmwrite([OUT.folder 'snow_density_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.snow_density]);
%             dlmwrite([OUT.folder 'ALT_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.ALT_max]);
%             dlmwrite([OUT.folder 'FT_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.FT_state]);
%             dlmwrite([OUT.folder 'snow_depth_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.snow_depth]);
%             dlmwrite([OUT.folder 'snow_covered_days_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.snow_covered_days]);
%             dlmwrite([OUT.folder 'T1m_max_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T1m_max]);
%             dlmwrite([OUT.folder 'T2m_max_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T2m_max]);
%             dlmwrite([OUT.folder 'T1m_min_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T1m_min]);
%             dlmwrite([OUT.folder 'T2m_min_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T2m_min]);
%             dlmwrite([OUT.folder 'T1m_std_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T1m_std]);
%             dlmwrite([OUT.folder 'T2m_std_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T2m_std]);
% 
% 	    stratigraphy = reshape(OUT.stratigraphy, size(OUT.stratigraphy,2)./OUT.ensemble_size,OUT.ensemble_size);
% 	    wind_compaction_timescale = reshape(OUT.wind_compaction_timescale, size(OUT.stratigraphy,2)./OUT.ensemble_size,OUT.ensemble_size);
% 	    sf_factor = reshape(OUT.sf_factor, size(OUT.stratigraphy,2)./OUT.ensemble_size,OUT.ensemble_size);
% 	    gain_loose = OUT.gain_loose;
%             save([OUT.folder 'stratigraphy_domain_' num2str(OUT.domain_id) '.mat'], 'stratigraphy', 'wind_compaction_timescale', 'sf_factor', 'gain_loose');
%             
%             if ~isempty(OUT.POI_id)
%                 dlmwrite([OUT.folder 'POI_domain_' num2str(OUT.domain_id) '.dat'], [OUT.POI.timestamp OUT.POI.GST OUT.POI.T1m OUT.POI.T2m OUT.POI.T10m OUT.POI.swe])
%             end
%                         
%             OUT.ACC.T1m=0;
%             OUT.ACC.T2m=0;
% 	    OUT.ACC.T5m=0;
% 	    OUT.ACC.T10m = 0;
% 	    OUT.ACC.T20m = 0;
% 
%             OUT.ACC.GST=0;
%             OUT.ACC.max_ALT = 0;
%             OUT.ACC.FT_state = 0;
%             OUT.ACC.N=0;
%             OUT.ACC.snow_depth = 0;
%             OUT.ACC.snow_covered_days = 0;
% 	    OUT.ACC.T1m_max = -1000;
% 	    OUT.ACC.T1m_min = 1000;
%             OUT.ACC.T2m_max = -1000;
% 	    OUT.ACC.T2m_min = 1000;
% 	    OUT.ACC.T1m_square = 0;
%             OUT.ACC.T2m_square = 0;
% 
%         end
%         
%     
% 
%         function OUT = get_yearly_av_OUT_forcing(OUT, forcing, profile, year)  %write averages of the forcing data
% 		
% 
% 	    OUT.surfT_av = [OUT.surfT_av; mean(forcing.surf_T,2)'];
%             OUT.surfT_FDD = [OUT.surfT_FDD; sum(-forcing.surf_T.*double(forcing.surf_T<0).*8, 2)'];
% 	    OUT.surfT_TDD = [OUT.surfT_TDD; sum(forcing.surf_T.*double(forcing.surf_T>0).*8, 2)'];
%             OUT.MODIS_weight = [OUT.MODIS_weight; mean(dlmread([forcing.folder 'forcing_MODIS_weight' num2str(year) '.dat']),2)']; 
% 	    OUT.theta_m = [OUT.theta_m;  mean(profile.theta_m(1:5,:),1)];
% 	    OUT.theta_o = [OUT.theta_o; mean(profile.theta_o(1:5,:),1)];
%             OUT.theta_w = [OUT.theta_w; mean(profile.theta_w(1:5,:),1)];
% 	    OUT.T_end_freezing = [OUT.T_end_freezing; mean(profile.T_end_freezing(1:5,:),1)];
%             
% 	    dlmwrite([OUT.folder 'surfT_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.surfT_av]);
% 	    dlmwrite([OUT.folder 'surf_FDD_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.surfT_FDD]);
% 	    dlmwrite([OUT.folder 'surf_TDD_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.surfT_TDD]);
% 	    dlmwrite([OUT.folder 'MODIS_weight_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.MODIS_weight]);
% 	    dlmwrite([OUT.folder 'theta_m_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.theta_m]);
% 	    dlmwrite([OUT.folder 'theta_o_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.theta_o]);
% 	    dlmwrite([OUT.folder 'theta_w_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.theta_w]);
%             dlmwrite([OUT.folder 'T_end_freezing_domain_' num2str(OUT.domain_id) '.dat'], [OUT.year_list OUT.T_end_freezing]);
% 
%         end
        
        
        function out = set2zero(out)
            
            out.TEMP.acc_T1m=single(0); %0; 
            out.TEMP.acc_T2m=single(0);
            out.TEMP.acc_T5m=single(0);
            out.TEMP.acc_T10m=single(0);
            out.TEMP.acc_T20m=single(0);
            out.TEMP.acc_GST=single(0);
            out.TEMP.acc_max_ALT = single(0);
            out.TEMP.acc_FT_state = 0;
            out.TEMP.acc_FT_isothermal = 1;
            out.TEMP.acc_N = 0;
            out.TEMP.acc_snow_density = single(0);
            out.TEMP.acc_snow_depth = single(0);
            out.TEMP.acc_snow_covered_days = single(0);
            out.TEMP.acc_T1m_max = single(-1000);
            out.TEMP.acc_T1m_min = single(1000);
            out.TEMP.acc_T2m_max = single(-1000);
            out.TEMP.acc_T2m_min = single(1000);
            out.TEMP.acc_T1m_square = single(0);
            out.TEMP.acc_T2m_square = single(0);
%             out.TEMP.acc_FDD = single(0);  %initialization
%             out.TEMP.acc_TDD = single(0);
        end
            
        function out = reset_STATVAR(out)
            out.STATVAR.T1m = [];
            out.STATVAR.T2m = [];
            out.STATVAR.T5m = [];
            out.STATVAR.T10m = [];
            out.STATVAR.T20m = [];
            out.STATVAR.GST = [];
            out.STATVAR.max_ALT  = [];
            out.STATVAR.FT_state = [];
            %out.STATVAR.FT_isothermal = [];
            out.STATVAR.snow_density = [];
            out.STATVAR.snow_depth = [];
            out.STATVAR.snow_covered_days  = [];
            out.STATVAR.T1m_max  = [];
            out.STATVAR.T1m_min = [];
            out.STATVAR.T2m_max = [];
            out.STATVAR.T2m_min = [];
            out.STATVAR.T1m_std = [];
            out.STATVAR.T2m_std = [];
%             out.STATVAR.FDD = []; %initialization
%             out.STATVAR.TDD = [];
            out.TEMP.scale_factor = [0.002;0.002;0.002; 0.002;0.002;0.002; 0.01; 1; 1; 0.01; 1; 0.002;0.002;0.002;0.002;0.002;0.002];
            out.TEMP.scale_offset = [-70;-70;-70;-70;-70;-70; 0; -10; 0; 0; 0; -70;-70;-70;-70;-70;-70];
        end

    end
end

