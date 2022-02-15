%========================================================================
% CryoGrid LATERAL_IA class LAT3D_SNOW_CROCUS 
% simulates lateral wind drift of snow between different CryoGrid
% stratigraphies. Blowing snow is assigned to stratigraphies with lower
% surface elevation.
% NOTE: works only with the SNOW classes SNOW_crocus_... and SNOW_crocus2_... 
% S. Westermann, Oct 2020
%========================================================================


classdef LAT3D_SNOW_CROCUS < BASE_LATERAL

    
    methods

        %----mandatory functions---------------
        %----initialization--------------------

        
        function lateral = provide_PARA(lateral)
            lateral.PARA.N_drift = []; %5;
            lateral.PARA.drift_loss_fraction = []; %0; %fraction of drifting snow that is lost from the system
            lateral.PARA.ia_time_increment = []; %0.05; %must be a multiple of the time increment of the main lateral class
%             lateral.PARA.ia_time_increment_min = []; %0.05;
%             lateral.PARA.ia_time_increment_max = []; %0.25;            
            %lateral.PARA.ia_time_next = [];
        end
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = []; %24 .* 3600;
        end
        
        
        function lateral = provide_STATVAR(lateral)
            
        end

        
        function lateral = finalize_init(lateral, tile)
            
        end
        
        %----time integration------------
        
        function lateral = pull(lateral, tile)
            lateral.PARENT.STATVAR2ALL.snow_drift = 0; % 0: no snow class; 2: driftable snow; 1: snow class, but snow is not driftable
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            while ~(strcmp(class(CURRENT), 'Bottom'))
                if strcmp(class(CURRENT), 'SNOW_crocus_bucketW_seb') || strcmp(class(CURRENT), 'SNOW_crocus2_bucketW_seb')
                    CURRENT = lateral3D_pull_snow(CURRENT, lateral);
                end
                CURRENT = CURRENT.NEXT;
            end
        end
        
        function lateral = get_derivatives(lateral, tile) %no need to loop through stratigraphy, all the information is in lateral.PARENT
            %calculate the exposure
            %loop over all ensemble members, mix the drifting part of the snow 
            
            lateral.STATVAR.snow_drift_yes_no = 0;
            if lateral.PARENT.STATVAR2ALL.snow_drift > 0 %uppermost class is SNOW
                altitude = lateral.PARENT.STATVAR2ALL.upperPos;
                area = lateral.PARENT.STATVAR2ALL.ds_area;
                lateral.STATVAR.snow_drift_yes_no = (lateral.PARENT.STATVAR2ALL.snow_drift == 2);
                %lateral.STATVAR.snow_drift_yes_no
                for j=1:size(lateral.PARENT.ENSEMBLE,1)
                    if lateral.PARENT.ENSEMBLE{j,1}.snow_drift >0
                        altitude = [altitude; lateral.PARENT.ENSEMBLE{j,1}.upperPos];
                        area = [area; lateral.PARENT.ENSEMBLE{j,1}.ds_area];
                        lateral.STATVAR.snow_drift_yes_no = lateral.STATVAR.snow_drift_yes_no + (lateral.PARENT.ENSEMBLE{j,1}.snow_drift == 2);
                    end
                end
                
                %lateral.STATVAR.snow_drift_yes_no
                if lateral.STATVAR.snow_drift_yes_no
                    
                    area_above = (altitude + 0.05 < altitude')*area;
                    area_below =(altitude - 0.05 > altitude')*area;
                    exposure2 = (area_above - area_below) ./ (area_above + area_below);
                    exposure2(isnan(exposure2)) = 0;
                    lateral.STATVAR.exposure = exposure2(1,1); %own exposure
                    exposure = exposure2(1,1);
                    k=0;
                    for j=1:size(lateral.PARENT.ENSEMBLE,1)
                        if lateral.PARENT.ENSEMBLE{j,1}.snow_drift >0
                            exposure = [exposure; exposure2(j+k+1,1)];
                        else
                            exposure = [exposure; 0];
                            k=k-1;
                        end
                    end
                    
                    
                    if lateral.STATVAR.exposure > 0 %own realization gains snow
                        area_acc = area(1,1) .* lateral.STATVAR.exposure;
                        lateral.STATVAR.ds.waterIce = 0;
                        lateral.STATVAR.ds.ice = 0;
                        lateral.STATVAR.ds.energy = 0;
                        %intensive variables - use waterIce as scaling variable,
                        %identical to ice when snow is driftable
                        lateral.STATVAR.ds.d = 0;
                        lateral.STATVAR.ds.s = 0;
                        lateral.STATVAR.ds.gs = 0;
                        lateral.STATVAR.ds.time_snowfall = 0;
                        lateral.STATVAR.ds.top_snow_date = 0;
                        lateral.STATVAR.ds.bottom_snow_date = 0;
                        volume=0;
                        for j=1:size(lateral.PARENT.ENSEMBLE,1)
                            
                                if lateral.PARENT.ENSEMBLE{j,1}.snow_drift > 1 && exposure(j+1,1) < 0 %all the loosing cells -> calculate total drifting snow pool
                                    volume = volume - exposure(j+1,1).*lateral.PARENT.ENSEMBLE{j,1}.ds_volume;
                                    lateral.STATVAR.ds.waterIce = lateral.STATVAR.ds.waterIce - exposure(j+1,1).* lateral.PARENT.ENSEMBLE{j,1}.ds_waterIce;
                                    lateral.STATVAR.ds.ice = lateral.STATVAR.ds.ice - exposure(j+1,1).* lateral.PARENT.ENSEMBLE{j,1}.ds_ice;
                                    lateral.STATVAR.ds.energy = lateral.STATVAR.ds.energy - exposure(j+1,1).* lateral.PARENT.ENSEMBLE{j,1}.ds_energy;
                                    %intensive variables - use waterIce as scaling variable,
                                    %identical to ice when snow is driftable
                                    lateral.STATVAR.ds.d = lateral.STATVAR.ds.d  - exposure(j+1,1).* lateral.PARENT.ENSEMBLE{j,1}.ds_waterIce .* lateral.PARENT.ENSEMBLE{j,1}.ds_d;
                                    lateral.STATVAR.ds.s = lateral.STATVAR.ds.s - exposure(j+1,1).* lateral.PARENT.ENSEMBLE{j,1}.ds_waterIce .*  lateral.PARENT.ENSEMBLE{j,1}.ds_s;
                                    lateral.STATVAR.ds.gs = lateral.STATVAR.ds.gs - exposure(j+1,1).* lateral.PARENT.ENSEMBLE{j,1}.ds_waterIce .* lateral.PARENT.ENSEMBLE{j,1}.ds_gs;
                                    lateral.STATVAR.ds.time_snowfall = lateral.STATVAR.ds.time_snowfall - exposure(j+1,1).* lateral.PARENT.ENSEMBLE{j,1}.ds_waterIce .* lateral.PARENT.ENSEMBLE{j,1}.ds_time_snowfall;
                                    
                                    lateral.STATVAR.ds.top_snow_date = lateral.STATVAR.ds.top_snow_date - exposure(j+1,1).* lateral.PARENT.ENSEMBLE{j,1}.ds_waterIce .* lateral.PARENT.ENSEMBLE{j,1}.ds_top_snow_date;
                                    lateral.STATVAR.ds.bottom_snow_date = lateral.STATVAR.ds.bottom_snow_date - exposure(j+1,1).* lateral.PARENT.ENSEMBLE{j,1}.ds_waterIce .* lateral.PARENT.ENSEMBLE{j,1}.ds_bottom_snow_date;
                                    
                                elseif lateral.PARENT.ENSEMBLE{j,1}.snow_drift >0 && exposure(j+1,1) > 0  %all the gaining cells ->
                                    area_acc = area_acc + area(j+1,1) .* exposure(j+1,1);
                                end
                            
                        end
                        lateral.STATVAR.ds.d = lateral.STATVAR.ds.d ./ lateral.STATVAR.ds.waterIce;
                        lateral.STATVAR.ds.d(isnan(lateral.STATVAR.ds.d)) = 0;
                        lateral.STATVAR.ds.s = lateral.STATVAR.ds.s ./ lateral.STATVAR.ds.waterIce;
                        lateral.STATVAR.ds.s(isnan(lateral.STATVAR.ds.s)) = 0;
                        lateral.STATVAR.ds.gs = lateral.STATVAR.ds.gs ./ lateral.STATVAR.ds.waterIce;
                        lateral.STATVAR.ds.gs(isnan(lateral.STATVAR.ds.gs)) = 0;
                        lateral.STATVAR.ds.time_snowfall = lateral.STATVAR.ds.time_snowfall ./ lateral.STATVAR.ds.waterIce;
                        lateral.STATVAR.ds.time_snowfall(isnan(lateral.STATVAR.ds.time_snowfall)) = 0;
                        
                        lateral.STATVAR.ds.top_snow_date = lateral.STATVAR.ds.top_snow_date./ lateral.STATVAR.ds.waterIce;
                        lateral.STATVAR.ds.top_snow_date(isnan(lateral.STATVAR.ds.top_snow_date)) = 0;
                        lateral.STATVAR.ds.bottom_snow_date = lateral.STATVAR.ds.bottom_snow_date./ lateral.STATVAR.ds.waterIce;
                        lateral.STATVAR.ds.bottom_snow_date(isnan(lateral.STATVAR.ds.bottom_snow_date)) = 0;
                        
                        gain_fraction = area(1,1) .* exposure(1,1) ./ area_acc;
                        volume = volume .* gain_fraction .* (1-lateral.PARA.drift_loss_fraction);
                        lateral.STATVAR.ds.waterIce = lateral.STATVAR.ds.waterIce .* gain_fraction .* (1-lateral.PARA.drift_loss_fraction);
                        lateral.STATVAR.ds.ice = lateral.STATVAR.ds.ice .* gain_fraction .* (1-lateral.PARA.drift_loss_fraction);
                        lateral.STATVAR.ds.energy = lateral.STATVAR.ds.energy .* gain_fraction .* (1-lateral.PARA.drift_loss_fraction);
                        lateral.STATVAR.ds.layerThick =  volume ./ area(1,1);
                        lateral.STATVAR.ds.target_density = lateral.STATVAR.ds.ice ./ volume;
                        lateral.STATVAR.ds.water = 0;
                    end
                end
            end
        end

        
        function lateral = push(lateral, tile)
            
            CURRENT = lateral.PARENT.TOP.NEXT; %find correct stratigraphy class
            while ~(strcmp(class(CURRENT), 'Bottom'))
                if strcmp(class(CURRENT), 'SNOW_crocus_bucketW_seb') || strcmp(class(CURRENT), 'SNOW_crocus2_bucketW_seb')
                    CURRENT = lateral3D_push_snow(CURRENT, lateral);
                end
                CURRENT = CURRENT.NEXT;
            end
        end
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 0;
            if t + lateral.PARENT.IA_TIME_INCREMENT >= lateral.PARA.ia_time_next - 1e-7
                lateral.PARENT.ACTIVE(i,1) = 1;
                %lateral.PARA.ia_time_next = t + lateral.PARENT.IA_TIME_INCREMENT + lateral.PARA.ia_time_increment;
                lateral.PARA.ia_time_next = lateral.PARA.ia_time_next + lateral.PARA.ia_time_increment;
            end
        end
        
        function lateral = set_ia_time(lateral, t)
            lateral.PARA.ia_time_next = t;
        end
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE_LATERAL(ground);
            
            ground.PARA.class_category = 'LATERAL_IA';
            
            ground.PARA.options = [];
            ground.PARA.STATVAR = [];

            ground.PARA.default_value.N_drift = {5};
            ground.PARA.comment.N_drift ={'factor translating CROCUS driftability index to snow removed per unit time'};
            
            ground.PARA.default_value.drift_loss_fraction = {0};
            ground.PARA.comment.drift_loss_fraction ={'fraction of drifting snow that is lost from the system'};
            
            ground.PARA.default_value.ia_time_increment = {0.25};
            ground.PARA.comment.ia_time_increment ={'time step [days], must be multiple of of LATERAL class timestep'};
        end
    end
    
end


