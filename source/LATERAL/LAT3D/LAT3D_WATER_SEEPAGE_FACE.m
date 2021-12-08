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
        
%         function lateral = LAT3D_WATER_SEEPAGE_FACE(index, pprovider, cprovider)
%             lateral@BASE_LATERAL(index, pprovider, cprovider);
%         end
        
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
                
            end
        end
        
        
        function lateral = push(lateral, tile)
            %taken care of by the 3D water exchange
        end
        
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 0;
            if t + lateral.PARENT.IA_TIME_INCREMENT >= lateral.PARA.ia_time_next - 1e-9
                lateral.PARENT.ACTIVE(i,1) = 1;
                lateral.PARA.ia_time_next = t + lateral.PARENT.IA_TIME_INCREMENT + lateral.PARA.ia_time_increment;
                %disp(lateral.PARA.ia_time_next-floor(lateral.PARA.ia_time_next));
            end
        end
        
        function lateral = set_ia_time(lateral, t)
            lateral.PARA.ia_time_next = t;
        end
        
    end
    
end


