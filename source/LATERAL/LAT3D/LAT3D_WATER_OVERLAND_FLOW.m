%========================================================================
% CryoGrid LATERAL_IA class LAT3D_WATER_OVERLAND_FLOW 

% NOTE: must be called after LAT3D_WATER_UNCONFINED_AQUIFER_OVERLAND_FLOW
% in the lateral class list!
% S. Westermann, Oct 2020
%========================================================================


classdef LAT3D_WATER_OVERLAND_FLOW < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = []; %24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.gradient = [];
            lateral.PARA.GaMa_coefficient = []; %Gauckler-Manning coefficient, https://en.wikipedia.org/wiki/Manning_formula
            lateral.PARA.overland_flow_contact_length = []; %40;
            lateral.PARA.overflow_threshold_elevation = [];
            lateral.PARA.ia_time_increment = []; %0.25; %must be a multiple of the time increment of the main lateral class
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.surface_run_off = [];
        end

        function lateral = finalize_init(lateral, tile)
            lateral.STATVAR.surface_run_off = 0;
        end
        
        %-----time integration-------
        
        function lateral = pull(lateral, tile)
            %taken care of by the 3D water exchange
        end
        
        function lateral = get_derivatives(lateral, tile) %no need to loop through stratigraphy, all the information is in lateral.PARENT
            %could be necessary to  limit flow for numerical stability?
            
            if lateral.PARENT.STATVAR.water_depth >1e-6 && ~isempty(lateral.PARENT.STATVAR.T_water) && lateral.PARENT.STATVAR.depths(1,1)>lateral.PARA.overflow_threshold_elevation
                water_depth = min(lateral.PARENT.STATVAR.water_depth, max(0, lateral.PARENT.STATVAR.depths(1,1) - lateral.PARA.overflow_threshold_elevation));
                velocity = lateral.PARA.GaMa_coefficient .* real(water_depth.^(2/3) .* abs(lateral.PARA.gradient).^0.5);
                flow = -velocity .* water_depth .* lateral.PARA.overland_flow_contact_length; %negative, outflow only
                flow = max(flow, -lateral.PARENT.STATVAR.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec));
                
                flow_energy = flow .* lateral.PARENT.CONST.c_w .* lateral.PARENT.STATVAR.T_water(1,1);
                
                lateral.PARENT.STATVAR.water_flux(1,1) = lateral.PARENT.STATVAR.water_flux(1,1) + flow .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
                lateral.PARENT.STATVAR.water_flux_energy(1,1) = lateral.PARENT.STATVAR.water_flux_energy(1,1) + flow_energy .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;       
                lateral.STATVAR.surface_run_off = lateral.STATVAR.surface_run_off - flow .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
                
                lateral.STATVAR.flow = flow;
                lateral.STATVAR.max_flow = -lateral.PARENT.STATVAR.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec);
                lateral.STATVAR.water_depth = lateral.PARENT.STATVAR.water_depth;
            end
            
        end
        
        
        function lateral = push(lateral, tile)
            %taken care of by the 3D water exchange
        end
        
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 0;
            if t + lateral.PARENT.IA_TIME_INCREMENT >= lateral.PARA.ia_time_next - 1e-7
                lateral.PARENT.ACTIVE(i,1) = 1;
                %lateral.PARA.ia_time_next = t + lateral.PARENT.IA_TIME_INCREMENT + lateral.PARA.ia_time_increment;
                lateral.PARA.ia_time_next = lateral.PARA.ia_time_next + lateral.PARA.ia_time_increment;
                %disp(lateral.PARA.ia_time_next-floor(lateral.PARA.ia_time_next));
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
            
            ground.PARA.default_value.gradient = {'0.01'};
            ground.PARA.comment.gradient = {'gradient of surface [vertial m/ horizontal m]'};
            
            ground.PARA.default_value.GaMa_coefficient = {15};
            ground.PARA.comment.GaMa_coefficient = {'Gauckler-Manning coefficient, https://en.wikipedia.org/wiki/Manning_formula'};
            
            ground.PARA.default_value.overland_flow_contact_length = {1};
            ground.PARA.comment.overland_flow_contact_length = {'lateral contact length for overland flow = width of channel [m]'};
            
            ground.PARA.default_value.overflow_threshold_elevation = {0};
            ground.PARA.comment.overflow_threshold_elevation = {'threshold elevation, no overland flow when water level is below [m a.s.l.]'};
            
            ground.PARA.default_value.ia_time_increment = {0.25};
            ground.PARA.comment.ia_time_increment ={'time step [days]'};
        end
    end
    
end


