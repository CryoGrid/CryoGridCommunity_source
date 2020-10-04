
classdef LAT3D_WATER_RESERVOIR < BASE_LATERAL

    
    methods
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = 24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.reservoir_elevation = 1221.1;
            lateral.PARA.reservoir_temperature = []; %only active for Xice classes - if empty, water added at the temperature of the respective grid cell
            lateral.PARA.hardBottom_cutoff = 0.03; %hard bottom if saturated and water content below
            lateral.PARA.distance_reservoir = 10; 
            lateral.PARA.reservoir_contact_length = 1;
            lateral.PARA.ia_time_increment = 0.25; %must be a multiple of the time increment of the main lateral class
            lateral.PARA.ia_time_next = [];
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.subsurface_run_off = [];
        end
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function lateral = finalize_init(lateral)
            lateral.STATVAR.subsurface_run_off = 0;
        end

        
        function lateral = pull(lateral)
            %taken care of by the 3D water exchange
        end
        
        function lateral = get_derivatives(lateral)
            
            %
            if lateral.PARENT.STATVAR.water_available
                head = -lateral.PARENT.STATVAR.water_table_elevation + lateral.PARA.reservoir_elevation;
                contact_height = lateral.PARENT.STATVAR.depths(1:end-1,1) - lateral.PARENT.STATVAR.depths(2:end,1);
                reservoir_flux = lateral.PARA.reservoir_contact_length .* contact_height .* lateral.PARENT.STATVAR.hydraulicConductivity .* head ./ lateral.PARA.distance_reservoir; 
                if isempty(lateral.PARA.reservoir_temperature) 
                        reservoir_flux_energy = reservoir_flux .* lateral.PARENT.CONST.c_w .* lateral.PARENT.STATVAR.T_water;
                else
                    reservoir_flux_energy = reservoir_flux .* lateral.PARENT.CONST.c_w .* ...
                        (double(reservoir_flux>0) .* lateral.PARA.reservoir_temperature + double(reservoir_flux<0) .* lateral.PARENT.STATVAR.T_water);
                end

                lateral.PARENT.STATVAR.water_flux = lateral.PARENT.STATVAR.water_flux + reservoir_flux .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
                lateral.PARENT.STATVAR.water_flux_energy = lateral.PARENT.STATVAR.water_flux_energy + reservoir_flux_energy .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;       
                lateral.STATVAR.subsurface_run_off = lateral.STATVAR.subsurface_run_off - sum(reservoir_flux) .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
                
            end
        end
        
        
        function lateral = push(lateral, forcing)
            %taken care of by the 3D water exchange
        end
        
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 0;
            if t + lateral.PARENT.IA_TIME_INCREMENT >= lateral.PARA.ia_time_next - 1e-9
                lateral.PARENT.ACTIVE(i,1) = 1;
                lateral.PARA.ia_time_next = t + lateral.PARENT.IA_TIME_INCREMENT + lateral.PARA.ia_time_increment;
            end
        end
        
        function lateral = set_ia_time(lateral, t)
            lateral.PARA.ia_time_next = t;
        end
        
        
        
    end
    
end


