classdef LATERAL_snow < LATERAL_water

    methods
        
        function xls_out = write_excel(~)
            xls_out = {'LATERAL','index',NaN,NaN;'LATERAL_snow',1,NaN,NaN;NaN,NaN,NaN,NaN;'interaction_timestep',1,'[hr]','interval for interaction between parallel realizations';'exposure',1,'[m]',NaN;'area',100,'[m^2]',NaN;'delta',0.100000000000000,'[m]','min difference in surface altitudes for snow exchange';'LATERAL_END',NaN,NaN,NaN};
        end
        
        function lateral = provide_variables(lateral)
            lateral = provide_variables@LATERAL_water(lateral);
            lateral.PARA.exposure = [];
            lateral.PARA.area = [];
            lateral.PARA.delta = [];
            lateral.PARA.N_drift = [];
        end
        
        function lateral = initalize_from_file(lateral, section)
            lateral = initalize_from_file@LATERAL_water(lateral, section);
        end
        
        function [lateral, forcing] = complete_init_lateral(lateral, forcing)
            [lateral, forcing] = complete_init_lateral@LATERAL_water(lateral, forcing);
        end
        
        function [lateral, snow] = lateral_interaction(lateral,snow,t)
            if t == lateral.INTERACTION_TIME
                [lateral, snow] = lateral_interaction@LATERAL_water(lateral,snow,t);
                labBarrier
                % Initialize interaction
                if strcmp(class(snow),'SNOW_simple_seb_crocus') || strcmp(class(snow),'SNOW_crocus')
                    lateral.STATUS.snow(labindex) = 1;
                else
                    lateral.STATUS.snow(labindex) = 0;
                end
                
                lateral.TEMP.exposures(labindex) = lateral.PARA.exposure(labindex) + sum(snow.STATVAR.layerThick);
                
                
                for j = 1:numlabs
                    if j ~= labindex
                        labSend(lateral.STATUS.snow(labindex),j,1);
                        labSend(lateral.TEMP.exposures(labindex),j,2);
                    end
                end
                for j = 1:numlabs
                    if j ~= labindex
                        lateral.STATUS.snow(j)      = labReceive(j,1);
                        lateral.TEMP.exposures(j)   = labReceive(j,2);
                    end
                end
                
                labBarrier
                if sum(lateral.STATUS.snow) >= 2
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
                    
%                     lateral.TEMP.surfaceAltitudes(labindex) = snow.STATVAR.upperPos;
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
                            labSend(snow_out,j,3);
%                             labSend(snow_out.water,j,4);
%                             labSend(snow_out.waterIce,j,5);
%                             labSend(snow_out.layerThick,j,6);
%                             labSend(snow_out.energy,j,7);
%                             labSend(snow_out.d,j,8);
%                             labSend(snow_out.s,j,9);
%                             labSend(snow_out.gs,j,10);
%                             labSend(snow_out.time_snowfall,j,11);
                        end
                    end
                    for j = 1:numlabs
                        if j ~= labindex
                            snow_out_j = labReceive(j,3);
                            lateral.TEMP.ice(j)             = snow_out_j.ice;
                            lateral.TEMP.water(j)           = snow_out_j.water;
                            lateral.TEMP.waterIce(j)        = snow_out_j.waterIce;
                            lateral.TEMP.layerThick(j)      = snow_out_j.layerThick;
                            lateral.TEMP.energy(j)          = snow_out_j.energy;
                            lateral.TEMP.d(j)               = snow_out_j.d;
                            lateral.TEMP.s(j)               = snow_out_j.s;
                            lateral.TEMP.gs(j)              = snow_out_j.gs;
                            lateral.TEMP.time_snowfall(j)   = snow_out_j.time_snowfall;
%                             lateral.TEMP.water(j)           = labReceive(j,4);
%                             lateral.TEMP.waterIce(j)        = labReceive(j,5);
%                             lateral.TEMP.layerThick(j)      = labReceive(j,6);
%                             lateral.TEMP.energy(j)          = labReceive(j,6);
%                             lateral.TEMP.d(j)               = labReceive(j,7);
%                             lateral.TEMP.s(j)               = labReceive(j,8);
%                             lateral.TEMP.gs(j)              = labReceive(j,9);
%                             lateral.TEMP.time_snowfall(j)   = labReceive(j,10);
                        end
                    end
 
                    % Add snow if deposition occurs
                    if drift_index(labindex) > 0 && sum(lateral.TEMP.ice) > 0
                        snow_drifting = get_snow_mixed(lateral);
                        snow_in = get_snow_deposited(drift_index,snow_drifting,lateral);
                        snow = add_drifting_snow(snow,snow_in);
                    end
                end
                
                labBarrier
            end
        end
    end
end

