%========================================================================
% CryoGrid LATERAL_IA class LAT_WATER_RESERVOIR 
% LAT_WATER_RESERVOIR simulates lateral water flow between a static water
% reservoir at defined  elevation (absolute elevation, not relative to the surface!), 
% as well as contact lengths (i.e. width) and distance from the GROUND column.
% Water temperatures of the reservoir can have a defined constant
% temperature (only relevant for inflow), otherwise inflow at grid cell temperature is assumed. 

% S. Westermann, Oct 2020
%========================================================================

classdef LAT_WATER_RESERVOIR < BASE_LATERAL
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        

        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = [];
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.reservoir_elevation = [];
            lateral.PARA.reservoir_temperature = []; %only active for Xice classes - if empty, water added at the temperature of the respective grid cell
            lateral.PARA.hardBottom_cutoff = []; %hard bottom if saturated and water content below
            lateral.PARA.distance_reservoir = []; 
            lateral.PARA.reservoir_contact_length = [];
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.subsurface_run_off = [];
        end
        
        function lateral = finalize_init(lateral, tile)
            lateral.STATVAR.subsurface_run_off = 0;
        end
        
        %-----time integration-------------
        
        %only push function needed
        function lateral = push(lateral, tile)
            CURRENT = lateral.PARENT.TOP.NEXT;
            lateral.TEMP.open_system = 1; %start with open system
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral_push_water_reservoir(CURRENT, lateral);
                CURRENT = compute_diagnostic(CURRENT, tile);
                CURRENT = CURRENT.NEXT;
            end
        end
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 1;
        end
        
        function lateral = get_derivatives(lateral, tile)
            
        end
        
        function lateral = pull(lateral, tile)
            
        end
        
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE_LATERAL(ground);
            
            ground.PARA.class_category = 'LATERAL_IA';
            
            ground.PARA.options = [];
            ground.PARA.STATVAR = [];
            
            ground.PARA.default_value.reservoir_temperature = {''};
            ground.PARA.comment.reservoir_temperature = {'temperature water in reservoir [degreeC], only active for Xice classes - if empty, water added at the temperature of the respective grid cell'};
            
            ground.PARA.default_value.reservoir_elevation = {19};
            ground.PARA.comment.reservoir_elevation = {'elevation of water reservoir [meter a.s.l.]'};
            
            ground.PARA.default_value.hardBottom_cutoff = {0.03};
            ground.PARA.comment.hardBottom_cutoff = {'hard bottom  = no water flow if saturated and water content below [vol. water content, -]'};

            ground.PARA.default_value.distance_reservoir = {10};
            ground.PARA.comment.distance_reservoir ={'distance to water reservoir [m]'};
            
            ground.PARA.default_value.reservoir_contact_length = {1};
            ground.PARA.comment.reservoir_contact_length = {'lateral contact length to water reservoir [m]'};
        end
    end
    
end


