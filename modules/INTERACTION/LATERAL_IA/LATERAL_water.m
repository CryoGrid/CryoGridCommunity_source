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
            lateral.PARA.hydraulic_conductivity = [];
            lateral.PARA.area = [];
            lateral.PARA.distance = [];
            lateral.PARA.contact_length = [];
            lateral.PARA.ghost = [];
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
            forcing.PARA.altitude = forcing.PARA.altitude + lateral.PARA.relative_elevation(labindex);
            
            interaction_timestep = nan(1,numlabs);
            interaction_timestep(labindex) = lateral.PARA.interaction_timestep;
            labBarrier
            for j = 1:numlabs
                if j ~= labindex
                    labSend(lateral.PARA.interaction_timestep,j,101);
                end
            end
            for j = 1:numlabs
                if j ~= labindex
                    interaction_timestep(j) = labReceive(j,101);
                end
            end
            
            if std(interaction_timestep) ~= 0
                error('Interaction times not equal or undefined')
            end
            
            lateral.STATUS.water = zeros(1,numlabs + lateral.PARA.ghost);
            lateral.STATUS.snow = zeros(1,numlabs);
            lateral.TEMP.lastWaterChange = nan(1,3);
            
            disp(['Elevation: ' num2str(forcing.PARA.altitude)  ' Area: ' num2str(lateral.PARA.area(labindex))])
            labBarrier
        end
        
        function [lateral, top_class] = lateral_interaction(lateral,top_class,t)
            if t == lateral.INTERACTION_TIME
                % Update interaction time
                [YY,MM,DD,HH,~,~] = datevec(lateral.INTERACTION_TIME);
                HH = round(HH + lateral.PARA.interaction_timestep);
                lateral.INTERACTION_TIME = datenum(YY,MM,DD,HH,0,0);
                
                
                lateral = precondition_lateral_exchange(lateral,top_class);
                
                %                 if  strcmp(class(top_class),'SNOW_crocus') % only for testing
                %                     lateral.STATUS.water(labindex) = 1;
                %                 else
                %                     lateral.STATUS.water(labindex) = 0;
                %                 end
                %
                %                 for j = 1:numlabs
                %                     if j ~= labindex
                %                         labSend(lateral.STATUS.water(labindex),j,102);
                %                     end
                %                 end
                %                 for j = 1:numlabs
                %                     if j ~= labindex
                %                         lateral.STATUS.water(j) = labReceive(j,102);
                %                     end
                %                 end
                % End testing
                
                labBarrier
                if sum(lateral.STATUS.water) >= 2
                    CLASS = class(top_class);
                    realizations = numlabs + lateral.PARA.ghost;
                    
                    pot_water_fluxes = zeros(realizations);
                    
                    lateral.TEMP.frost_table  = zeros(1,realizations);
                    lateral.TEMP.water_table  = zeros(1,realizations);
                    lateral.TEMP.mobile_water = zeros(1,realizations);
%                     
                    if strcmp(CLASS(1:4),'SNOW')
                        [exchange, mobile_water] = get_water_exchange_SNOW(top_class);
                    elseif strcmp(CLASS(1:6),'GROUND')
                        [exchange, mobile_water] = get_water_exchange_GROUND(top_class);
                    end
                    
                    lateral.TEMP.frost_table(labindex)  = exchange.frost_table;
                    lateral.TEMP.water_table(labindex)  = exchange.water_table;
                    lateral.TEMP.mobile_water(labindex) = exchange.mobile_water;
                    
                    lateral.PARA.hydraulic_conductivity(labindex) = top_class.PARA.hydraulicConductivity;
                    
                    for j = 1:numlabs
                        if j ~= labindex && lateral.STATUS.water(j) == 1
                            labSend(exchange,j,104);
                            labSend(top_class.PARA.hydraulicConductivity,j,105);
                        end
                    end
                    for j = 1:numlabs
                        if j ~= labindex && lateral.STATUS.water(j) == 1
                            exchange_j = labReceive(j,104);
                            lateral.TEMP.frost_table(j)  = exchange_j.frost_table;
                            lateral.TEMP.water_table(j)  = exchange_j.water_table;
                            lateral.TEMP.mobile_water(j) = exchange_j.mobile_water;
                            lateral.PARA.hydraulic_conductivity(j) = labReceive(j,105);
                            
                            pot_water_fluxes = lateral_Darcy_flux(lateral,pot_water_fluxes,j);
                        end
                    end
                    % Boundary water flux
                    if lateral.PARA.ghost == 1
                        lateral.TEMP.frost_table(end) = 299;
                        lateral.TEMP.water_table(end) = 299;
                        lateral.TEMP.mobile_water(end) = 0;
                        pot_water_fluxes = lateral_Darcy_flux(lateral,pot_water_fluxes,numlabs+1);
                    end

                    water_fluxes_index = available_lateral_water(lateral,pot_water_fluxes);
                    water_fluxes_ensamble = zeros(realizations,realizations,realizations);
                    water_fluxes_ensamble(:,:,labindex) = water_fluxes_index;
                    
                    for j = 1:numlabs
                        if j ~= labindex && lateral.STATUS.water(j) == 1
                            labSend(water_fluxes_index,j,106);
                        end
                    end
                    for j = 1:numlabs
                        if j ~= labindex && lateral.STATUS.water(j) == 1
                            water_fluxes_j = labReceive(j,106);
                            water_fluxes_ensamble(:,:,j) = water_fluxes_j;
                        end
                    end
                    %
                    water_fluxes = nansum(water_fluxes_ensamble,3);
                    waterflux = nansum(water_fluxes(:,labindex)); % + boundary weater (implement later)
                    [lateral, waterflux] = lateral_water_oscillations(lateral, waterflux);
                    
                    if waterflux < 0    % loosing water
                        top_class = drain_water(top_class,mobile_water,waterflux);
                    elseif waterflux > 0 % gaining water
                        if strcmp(CLASS(1:4),'SNOW')
                            top_class = add_water_SNOW(top_class,waterflux);
                        elseif strcmp(CLASS(1:6),'GROUND')
                            top_class = add_water_GROUND(top_class,waterflux);
                        end
                    end
                    
                end
            end
        end
        
    end
end