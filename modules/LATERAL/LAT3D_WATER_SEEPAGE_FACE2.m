
classdef LAT3D_WATER_SEEPAGE_FACE2 < BASE_LATERAL

    
    methods
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = 24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.upperElevation = Inf;
            lateral.PARA.lowerElevation = 16;
            lateral.PARA.hardBottom_cutoff = 0.03; %hard bottom if saturated and water content below
            lateral.PARA.distance_seepageFace = 10; 
            lateral.PARA.seepage_contact_length = 1;
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
            number_of_own_aquifers = max(lateral.PARENT.STATVAR.aquifer_index);
            for j=1:number_of_own_aquifers
                a=find(lateral.PARENT.STATVAR.aquifer_index==j, 1, 'first');
                b=find(lateral.PARENT.STATVAR.aquifer_index==j, 1, 'last');
                
                cell_1 = lateral.PARENT.STATVAR.depths(a:b+1,1);
                cell_2 = [lateral.PARA.upperElevation; lateral.PARA.lowerElevation];
                overlap = get_overlap_aquifers(lateral.PARENT, cell_1, cell_2);
                overlap=[overlap zeros(size(overlap,1),3)];
                
                effective_hydraulic_conductivity=0;
                effective_hydraulic_conductivity_times_head = 0;
                effective_hydraulic_conductivity_times_head_reservoir = 0;
                for l=1:size(overlap,1)
                    ehc = lateral.PARENT.STATVAR.hydraulicConductivity(overlap(l,1)+a-1,1) ./ lateral.PARA.distance_seepageFace;
                    ehc = ehc .* lateral.PARA.seepage_contact_length;
                    overlap(l,5) = ehc;
                    overlap(l,6) = ehc .* (double(lateral.PARENT.STATVAR.water_status(overlap(l,1)+a-1,1)>=0) .* lateral.PARENT.STATVAR.head(overlap(l,1)+a-1,1) + ...
                        double(lateral.PARENT.STATVAR.water_status(overlap(l,1)+a-1,1)<0) .* overlap(l,4));
                    overlap(l,7) = ehc .* overlap(l,4);  %chnage this for ssepage face, rest should be fine!
                    
                    effective_hydraulic_conductivity =  effective_hydraulic_conductivity + ehc;
                    effective_hydraulic_conductivity_times_head = effective_hydraulic_conductivity_times_head + overlap(l,6);
                    effective_hydraulic_conductivity_times_head_reservoir = effective_hydraulic_conductivity_times_head_reservoir + overlap(l,7);
                end
                lateral.PARENT.STATVAR.external_fluxes = [lateral.PARENT.STATVAR.external_fluxes; ...
                    [j; effective_hydraulic_conductivity; effective_hydraulic_conductivity_times_head; effective_hydraulic_conductivity_times_head_reservoir]];
                
                lateral.PARENT.STATVAR_PRIVATE.overlap_info = [lateral.PARENT.STATVAR_PRIVATE.overlap_info; overlap];
                lateral.PARENT.STATVAR_PRIVATE.overlap_aquifer_index = [lateral.PARENT.STATVAR_PRIVATE.overlap_aquifer_index; j];
            end
        end
        
        function lateral = get_derivatives(lateral) %no need to loop through stratigraphy, all the information is in lateral.PARENT
            
            %seepage flow of saturated cells above other cell's water table
%             if lateral.PARENT.STATVAR.water_available
%                 depth_rel2surface = lateral.PARENT.STATVAR.water_table_elevation - lateral.PARENT.STATVAR.depths;
%                 depth_rel2surface = max(0,depth_rel2surface);
%                 head = (depth_rel2surface(1:end-1,1) + depth_rel2surface(2:end,1))/2;
%                 
%                 depths = min(max(lateral.PARENT.STATVAR.depths, lateral.PARA.lowerElevation), lateral.PARA.upperElevation);
%                 contact_height = (depths(1:end-1,1) - depths(2:end,1));
%                 
%                 seepage_flux =  - lateral.PARA.seepage_contact_length .* contact_height .* lateral.PARENT.STATVAR.hydraulicConductivity .* head ./ lateral.PARA.distance_seepageFace; %outflow!
%                 seepage_flux_energy = seepage_flux .* lateral.PARENT.CONST.c_w .* lateral.PARENT.STATVAR.T_water;
% 
%                 lateral.PARENT.STATVAR.water_flux = lateral.PARENT.STATVAR.water_flux + seepage_flux .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
%                 lateral.PARENT.STATVAR.water_flux_energy = lateral.PARENT.STATVAR.water_flux_energy + seepage_flux_energy .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;       
%                 lateral.STATVAR.subsurface_run_off = lateral.STATVAR.subsurface_run_off - sum(seepage_flux) .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
%                 
%             end
        end
        
        
        function lateral = push(lateral, forcing)
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


