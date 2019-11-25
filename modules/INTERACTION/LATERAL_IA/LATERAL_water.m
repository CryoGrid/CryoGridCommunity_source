classdef LATERAL_water
    
    properties
        PARA
        TEMP
        STATUS
        INTERACTION_TIME
    end
    
    methods
        
        function xls_out = write_excel(out)
            xls_out = {'LATERAL','index',NaN,NaN;'LATERAL_water',1,NaN,NaN;NaN,NaN,NaN,NaN;'interaction_timestep',1,'[hr]','interval for interaction between parallel realizations';'exposure',1,'[m]',NaN;'area',100,'[m^2]',NaN;'delta',0.100000000000000,'[m]','min difference in surface altitudes for snow exchange';'LATERAL_END',NaN,NaN,NaN};
        end
        
        function lateral = provide_variables(lateral)
            lateral.PARA.interaction_timestep = [];
            lateral.PARA.relative_elevation = [];
            lateral.PARA.area = [];
            lateral.PARA.distance = [];
            lateral.PARA.contact_length = [];
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
            lateral.STATUS.water = zeros(1,numlabs);
            
            disp(['Elevation: ' num2str(forcing.PARA.altitude)  ' Area: ' num2str(area(labindex))])
            labBarrier
        end
        
        function [lateral, top_class] = lateral_interaction(lateral,top_class,t)
            if t == lateral.INTERACTION_TIME
                % Update interaction time
                [YY,MM,DD,HH,~,~] = datevec(lateral.INTERACTION_TIME);
                HH = round(HH + lateral.PARA.interaction_timestep);
                lateral.INTERACTION_TIME = datenum(YY,MM,DD,HH,0,0);
                
%                 for j = 1:numlabs
%                     if j ~= labindex
%                         labSend(class(top_class),j,102);
%                     end
%                 end
%                 for j = 1:numlabs
%                     if j ~= labindex
%                         class_j = labReceive(j,102);
%                         if strcmp(class_j,class(top_class))
%                             lateral.STATUS.water(j) = 1;
%                             lateral.STATUS.water(labindex) = 1;
%                         end
%                     end
%                 end
                % For testing:
                if  strcmp(class(top_class),'SNOW_simple_seb_crocus') % only for testing
                    lateral.STATUS.water(labindex) = 1;
                else
                    lateral.STATUS.water(labindex) = 0;
                end
                for j = 1:numlabs
                    if j ~= labindex
                        labSend(lateral.STATUS.water(labindex),j,102);
                    end
                end
                for j = 1:numlabs
                    if j ~= labindex
                        lateral.STATUS.water(j) = labReceive(j,102);
                    end
                end
                % End testing
                
                labBarrier
                if sum(lateral.STATUS.water) >= 2
                    pot_water_fluxes = zeros(numlabs);

                    lateral.TEMP.frost_table  = zeros(1,numlabs);
                    lateral.TEMP.water_table  = zeros(1,numlabs);
                    lateral.TEMP.mobile_water = zeros(1,numlabs);
                    
                    % Adapt mobile water for GROUND later
                    mobile_water = top_class.STATVAR.water - (top_class.STATVAR.layerThick - top_class.STATVAR.ice).*top_class.PARA.field_capacity;
                    i = find(mobile_water > .001,1,'first'); % min threshld of mobile water for lateral fluxes
                    j = find(mobile_water > .001,1,'last'); % last cell with mobile water
                    mobile_water(mobile_water <= .001) = 0;
                    if sum(mobile_water) > 0
                        exchange.water_table = top_class.STATVAR.lowerPos + sum(top_class.STATVAR.layerThick(i+1:end)) + mobile_water(i)/(top_class.STATVAR.layerThick(i) - top_class.STATVAR.ice(i))*top_class.STATVAR.layerThick(i);
                    else
                        exchange.water_table = top_class.STATVAR.lowerPos;
                    end
                    exchange.frost_table = top_class.STATVAR.lowerPos + sum(top_class.STATVAR.layerThick(j+1:end));
                    exchange.mobile_water = sum(mobile_water);
                    
                    lateral.TEMP.frost_table(labindex)  = exchange.frost_table;
                    lateral.TEMP.water_table(labindex)  = exchange.water_table;
                    lateral.TEMP.mobile_water(labindex) = exchange.mobile_water;
                    
                    
                    for j = 1:numlabs
                        if j ~= labindex
                            labSend(exchange,j,103);
                        end
                    end
                    for j = 1:numlabs
                        if j ~= labindex
                            exchange_j = labReceive(j,103);
                            lateral.TEMP.frost_table(j)  = exchange_j.frost_table;
                            lateral.TEMP.water_table(j)  = exchange_j.water_table;
                            lateral.TEMP.mobile_water(j) = exchange_j.mobile_water;
                            
                            pot_water_fluxes = lateral_Darcy_flux(lateral,pot_water_fluxes,j);
                        end
                    end
                    pot_water_fluxes = lateral_Darcy_flux(lateral,pot_water_fluxes,2);
                    water_fluxes_index = available_lateral_water(lateral,pot_water_fluxes);
                    water_fluxes_ensamble = zeros(numlabs,numlabs,numlabs);
                    water_fluxes_ensamble(:,:,labindex) = water_fluxes_index;
                    
                    for j = 1:numlabs
                        if j ~= labindex
                            labSend(water_fluxes_index,j,106);
                        end
                    end
                    for j = 1:numlabs
                        if j ~= labindex
                            water_fluxes_j = labReceive(j,106);
                            water_fluxes_ensamble(:,:,j) = water_fluxes_j;
                        end
                    end
%                     
                    water_fluxes = nansum(water_fluxes_ensamble,3);
                    waterflux = nansum(water_fluxes(:,labindex)); % + boundary weater (implement later)
                    % only for testing
                    
                    if waterflux < 0    % loosing water
                        top_class = drain_water(top_class,mobile_water,waterflux);
                    elseif waterflux > 0 % gaining water
                        top_class = add_water(top_class,waterflux);
                    end
                    
                end
            end
        end
        
    end
end