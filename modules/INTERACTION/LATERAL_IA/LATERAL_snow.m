classdef LATERAL_snow
    
    properties
        PARA
        TEMP
        STATUS
        INTERACTION_TIME
    end
    
    methods
        
        function xls_out = write_excel(out)
            xls_out = {'LATERAL','index',NaN,NaN;'LATERAL_snow',1,NaN,NaN;NaN,NaN,NaN,NaN;'interaction_timestep',1,'[hr]','interval for interaction between parallel realizations';'relative_elevation',1,'[m]',NaN;'area',100,'[m^2]',NaN;'delta',0.100000000000000,'[m]','min difference in surface altitudes for snow exchange';'LATERAL_END',NaN,NaN,NaN};
        end
        
        function lateral = provide_variables(lateral)
            lateral.PARA.interaction_timestep = [];
            lateral.PARA.relative_elevation = [];
            lateral.PARA.area = [];
            lateral.PARA.delta = [];
        end
        
        function lateral = initalize_from_file(lateral, section)
            variables = fieldnames(lateral.PARA);
            for i=1:size(variables,1)
                for j=1:size(section,1)
                    if strcmp(variables{i,1}, section{j,1})
                        lateral.PARA.(variables{i,1}) = section{j,2};
                    end
                end
            end
        end
        
        function [lateral, forcing] = complete_init_lateral(lateral, forcing)
            lateral.INTERACTION_TIME = forcing.PARA.start_time + lateral.PARA.interaction_timestep/24;
            forcing.PARA.altitude = forcing.PARA.altitude + lateral.PARA.relative_elevation;
            
            area = nan(1,numlabs);
            area(labindex) = lateral.PARA.area;
            interaction_timestep = nan(1,numlabs);
            interaction_timestep(labindex) = lateral.PARA.interaction_timestep;
            labBarrier
            for j = 1:numlabs
                if j ~= labindex
                    labSend(lateral.PARA.area,j,100);
                    labSend(lateral.PARA.interaction_timestep,j,101);
                end
            end
            for j = 1:numlabs
                if j ~= labindex
                    area(j) = labReceive(j,100);
                    interaction_timestep(j) = labReceive(j,101);
                end
            end
            
            if std(interaction_timestep) ~= 0
                error('Interaction times not equal or undefined')
            end
            
            lateral.PARA.area = area;
            lateral.STATUS = zeros(1,numlabs);
            
            disp(['Elevation: ' num2str(forcing.PARA.altitude) ' Area: ' num2str(area(labindex))])
            labBarrier
        end
        
        function [lateral, snow] = lateral_interaction(lateral,snow,t)
            if t == lateral.INTERACTION_TIME
                % Update interaction time
                [YY,MM,DD,HH,~,~] = datevec(lateral.INTERACTION_TIME);
                HH = round(HH + lateral.PARA.interaction_timestep);
                lateral.INTERACTION_TIME = datenum(YY,MM,DD,HH,0,0);
                
                % Initialize interaction
                if strcmp(class(snow),'SNOW_simple_seb_crocus')
                    lateral.STATUS(labindex) = 1;
                else
                    lateral.STATUS(labindex) = 0;
                end
                
                lateral.TEMP.surfaceAltitudes(labindex) = snow.STATVAR.upperPos;
                
                
                for j = 1:numlabs
                    if j ~= labindex
                        labSend(lateral.STATUS(labindex),j,102);
                        labSend(snow.STATVAR.upperPos,j,103);
                    end
                end
                for j = 1:numlabs
                    if j ~= labindex
                        lateral.STATUS(j)                   = labReceive(j,102);
                        lateral.TEMP.surfaceAltitudes(j)    = labReceive(j,103);
                    end
                end
                
                labBarrier
                if sum(lateral.STATUS) >= 2
                    % determine exchange coefficient
                    drift_index = drift_exchange_index(lateral);
                    
                    
                    % Remove snow if erosion occurs
                    if drift_index(labindex) < 0
                        [snow, snow_out] = get_snow_eroded2(snow,lateral,-drift_index(labindex));
                        while snow.STATVAR.ice(1) == 0 % Whole uppermost layer was eroded
                            snow.STATVAR.ice(1)         = [];
                            snow.STATVAR.water(1)       = [];
                            snow.STATVAR.waterIce(1)    = [];
                            snow.STATVAR.layerThick(1)  = [];
                            snow.STATVAR.energy(1)      = [];
                            snow.STATVAR.d(1)           = [];
                            snow.STATVAR.s(1)           = [];
                            snow.STATVAR.gs(1)          = [];
                            snow.STATVAR.time_snowfall(1)   = [];
                            snow.STATVAR.target_density(1)  = [];
                        end
                    else
                        snow_out.ice            = 0;
                        snow_out.water          = 0;
                        snow_out.waterIce       = 0;
                        snow_out.layerThick     = 0;
                        snow_out.energy         = 0;
                        snow_out.d              = 0;
                        snow_out.s              = 0;
                        snow_out.gs             = 0;
                        snow_out.time_snowfall  = 0;
                    end
                    
                    lateral.TEMP.surfaceAltitudes(labindex) = snow.STATVAR.upperPos;
                    lateral.TEMP.ice(labindex)              = snow_out.ice;
                    lateral.TEMP.water(labindex)            = snow_out.water;
                    lateral.TEMP.waterIce(labindex)         = snow_out.waterIce;
                    lateral.TEMP.layerThick(labindex)       = snow_out.layerThick;
                    lateral.TEMP.energy(labindex)           = snow_out.energy;
                    lateral.TEMP.d(labindex)                = snow_out.d;
                    lateral.TEMP.s(labindex)                = snow_out.s;
                    lateral.TEMP.gs(labindex)               = snow_out.gs;
                    lateral.TEMP.time_snowfall(labindex)    = snow_out.time_snowfall;
                    
                    % Exchange snow properties
                    for j = 1:numlabs
                        if j ~= labindex
                            labSend(snow_out.ice,j,2);
                            labSend(snow_out.water,j,3);
                            labSend(snow_out.waterIce,j,4);
                            labSend(snow_out.layerThick,j,5);
                            labSend(snow_out.energy,j,6);
                            labSend(snow_out.d,j,7);
                            labSend(snow_out.s,j,8);
                            labSend(snow_out.gs,j,9);
                            labSend(snow_out.time_snowfall,j,10);
                        end
                    end
                    for j = 1:numlabs
                        if j ~= labindex
                            lateral.TEMP.ice(j)             = labReceive(j,2);
                            lateral.TEMP.water(j)           = labReceive(j,3);
                            lateral.TEMP.waterIce(j)        = labReceive(j,4);
                            lateral.TEMP.layerThick(j)      = labReceive(j,5);
                            lateral.TEMP.energy(j)          = labReceive(j,6);
                            lateral.TEMP.d(j)               = labReceive(j,7);
                            lateral.TEMP.s(j)               = labReceive(j,8);
                            lateral.TEMP.gs(j)              = labReceive(j,9);
                            lateral.TEMP.time_snowfall(j)   = labReceive(j,10);
                        end
                    end
                    
                    % Add snow if deposition occurs
                    if drift_index(labindex) > 0 && sum(lateral.TEMP.ice) > 0
                        snow_drifting = get_snow_mixed(lateral);
                        snow_in = get_snow_deposited(drift_index,snow_drifting);
                        snow = add_drifting_snow(snow,snow_in);
                    end
                end
                
                labBarrier
            end
        end
    end
end

