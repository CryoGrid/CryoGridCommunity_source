%========================================================================
% CryoGrid LATERAL_IA class LAT3D_WATER_SEEPAGE_FACE 
% simulates lateral water flow through a seepage face with defined upper
% and lower elevation (absolute elevation, not relative to the surface!), 
%  as well as contact lengths (i.e. width) and distance from the GROUND column.
% At the seepage face, air pressure is assumed.
% NOTE: must be called after LAT3D_WATER_UNCONFINED_AQUIFER or LAT3D_WATER
% in the lateral class list!
% S. Westermann, Oct 2020
%========================================================================


classdef LAT3D_WATER_SEEPAGE_FACE < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = []; %24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.upperElevation = []; %Inf;
            lateral.PARA.lowerElevation = []; %10;
            lateral.PARA.hardBottom_cutoff = []; %0.03; %hard bottom if saturated and water content below
            lateral.PARA.distance_seepageFace = []; %10; 
            lateral.PARA.seepage_contact_length = []; %40;
            lateral.PARA.ia_time_increment = []; %0.25; %must be a multiple of the time increment of the main lateral class
            %lateral.PARA.ia_time_next = [];
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.subsurface_run_off = [];
        end

        function lateral = finalize_init(lateral, tile)
            lateral.STATVAR.subsurface_run_off = 0;
        end
        
        %-----time integration-------
        
        function lateral = pull(lateral, tile)
            %taken care of by the 3D water exchange
        end
        
        function lateral = get_derivatives(lateral, tile) %no need to loop through stratigraphy, all the information is in lateral.PARENT
            
            %disp(['Hallo ' datestr(tile.t)])
            
            %seepage flow of saturated cells above other cell's water table
            if lateral.PARENT.STATVAR.water_available
                depth_rel2surface = lateral.PARENT.STATVAR.water_table_elevation - lateral.PARENT.STATVAR.depths;
                depth_rel2surface = max(0,depth_rel2surface);
                head = (depth_rel2surface(1:end-1,1) + depth_rel2surface(2:end,1))/2;
                
                depths = min(max(lateral.PARENT.STATVAR.depths, lateral.PARA.lowerElevation), lateral.PARA.upperElevation);
                contact_height = (depths(1:end-1,1) - depths(2:end,1));
                
                seepage_flux =  - lateral.PARA.seepage_contact_length .* contact_height .* lateral.PARENT.STATVAR.hydraulicConductivity .* head ./ lateral.PARA.distance_seepageFace; %outflow!
                seepage_flux_energy = seepage_flux .* lateral.PARENT.CONST.c_w .* lateral.PARENT.STATVAR.T_water;

                lateral.PARENT.STATVAR.water_flux = lateral.PARENT.STATVAR.water_flux + seepage_flux .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
                lateral.PARENT.STATVAR.water_flux_energy = lateral.PARENT.STATVAR.water_flux_energy + seepage_flux_energy .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;       
                lateral.STATVAR.subsurface_run_off = lateral.STATVAR.subsurface_run_off - sum(seepage_flux) .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
                
                %remove later

                lateral.STATVAR.head_test = head;
                lateral.STATVAR.depths = depths;
                lateral.STATVAR.contact_height = contact_height;
                lateral.STATVAR.seepage_flux = seepage_flux;

            end
            lateral.STATVAR.hydCond = lateral.PARENT.STATVAR.hydraulicConductivity;
            lateral.STATVAR.active = lateral.PARENT.ACTIVE;
            lateral.STATVAR.wt_elevation = lateral.PARENT.STATVAR.water_table_elevation;
            lateral.STATVAR.w_avail = lateral.PARENT.STATVAR.water_available;
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
            
            ground.PARA.default_value.upperElevation = {10000};
            ground.PARA.comment.upperElevation = {'upper elevation of seepage face [meter a.s.l.]'};
            
            ground.PARA.default_value.lowerElevation = {0};
            ground.PARA.comment.lowerElevation = {'lower elevation of seepage face [meter a.s.l.]'};
            
            ground.PARA.default_value.hardBottom_cutoff = {0.03};
            ground.PARA.comment.hardBottom_cutoff = {'hard bottom  = no water flow if saturated and water content below [vol. water content, -]'};
            
            ground.PARA.default_value.distance_seepageFace = {1};
            ground.PARA.comment.distance_seepageFace ={'distance to seepage face [m]'};
            
            ground.PARA.default_value.seepage_contact_length = {1};
            ground.PARA.comment.seepage_contact_length = {'lateral contact length (=width) of seepage face [m]'};
            
            ground.PARA.default_value.ia_time_increment = {0.25};
            ground.PARA.comment.ia_time_increment ={'time step [days], must be multiple of LATERAL class timestep'};
        end
        
    end
    
end


