%========================================================================
% CryoGrid LATERAL_IA class LAT3D_WATER_UNCONFINED_AQUIFER 
% simulates lateral water flow between pairs of CryoGrid stratigraphies
% for the topmost unconfined aquifer.
% NOTE: no flow if there is no unconfined aquifer for one of the two stratigraphies,
% e.g. if the first cell is saturated with ice. Use LAT_3D_WATER instead. 
% S. Westermann, Oct 2020
%========================================================================


classdef LAT3D_WATER_UNCONFINED_AQUIFER2 < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = []; %24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.hardBottom_cutoff = []; %0.03; %hard bottom if saturated and water content below
            %lateral.PARA.distance_reservoir = []; %2; 
           % lateral.PARA.reservoir_contact_length = []; %1;
            lateral.PARA.ia_time_increment = []; %0.25; %must be a multiple of the time increment of the main lateral class
            %lateral.PARA.ia_time_next = [];
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.subsurface_run_off = [];
            lateral.STATVAR.surface_runoff = [];
        end

        function lateral = finalize_init(lateral, tile)
            lateral.STATVAR.subsurface_run_off = 0;
            lateral.STATVAR.surface_runoff = 0;
        end
        
        %--- time integration---------------
        
        function lateral = pull(lateral, tile)
            
            lateral.PARENT.STATVAR.depths = [];
            lateral.PARENT.STATVAR.water_status = [];
            lateral.PARENT.STATVAR.hydraulicConductivity = [];
            lateral.PARENT.STATVAR.water_table_elevation = [];
            lateral.PARENT.STATVAR.water_available = 0;
            lateral.PARENT.STATVAR.T_water = [];
            lateral.PARENT.STATVAR.water_table_top_cell = 0;         
            
            lateral.TEMP.open_system = 1; %start with open system

            CURRENT = lateral.PARENT.TOP.NEXT;
            while ~(strcmp(class(CURRENT), 'Bottom')) && lateral.TEMP.open_system == 1
                CURRENT = lateral3D_pull_water_unconfined_aquifer(CURRENT, lateral);
                CURRENT = CURRENT.NEXT;
            end
        end
        
        function lateral = get_derivatives(lateral, tile) %no need to loop through stratigraphy, all the information is in lateral.PARENT

            lateral.PARENT = get_overlap_cells2(lateral.PARENT, 'depths', 'overlap');
            
            %calculate fluxes
            flux = lateral.PARENT.STATVAR.hydraulicConductivity .* 0;
            flux_energy = flux;
            for j=1:size(lateral.PARENT.ENSEMBLE,1)
                if lateral.PARENT.PARA.connected(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index)  %add if water is available at all
                    %flow between saturated cells
                    contact_length = lateral.PARENT.PARA.contact_length(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
                    distance = lateral.PARENT.PARA.distance(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
                    for i=1:size(lateral.PARENT.ENSEMBLE{j,1}.overlap,1) %does not do anything when overlap is empty!
                        
                        hc1 = lateral.PARENT.STATVAR.hydraulicConductivity(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1);
                        hc2 = lateral.PARENT.ENSEMBLE{j,1}.hydraulicConductivity(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,2),1);
                        hydraulicConductivity = hc1 .* hc2 ./ (hc1 .* distance./2 + hc2 .* distance./2); %change to different distances laters
                        hydraulicConductivity(isnan(hydraulicConductivity)) = 0;

                        flux_i = contact_length .* lateral.PARENT.ENSEMBLE{j,1}.overlap(i,3) .* hydraulicConductivity .* ...
                            -(lateral.PARENT.STATVAR.water_table_elevation - lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation);
                        flux(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) = flux(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) + flux_i;
                        flux_energy(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) = flux_energy(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) + flux_i .* lateral.PARENT.CONST.c_w .* ...
                            (double(flux_i<0).*lateral.PARENT.STATVAR.T_water(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) + ...
                            double(flux_i>=0) .* lateral.PARENT.ENSEMBLE{j,1}.T_water(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,2),1));
                        
                    end
                    %seepage flow of saturated cells above other cell's water table
                    if lateral.PARENT.STATVAR.water_available || lateral.PARENT.ENSEMBLE{j,1}.water_available
                        if lateral.PARENT.STATVAR.water_table_elevation < lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation %inflow
                            depth_rel2surface = lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation - lateral.PARENT.ENSEMBLE{j,1}.depths;
                            depth2other_worker_water_table = lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation - lateral.PARENT.STATVAR.water_table_elevation;
                            depth_rel2surface = max(0,depth_rel2surface);
                            depth_rel2surface = min(depth2other_worker_water_table,depth_rel2surface);
                            %head = (depth_rel2surface(1:end-1,1) + depth_rel2surface(2:end,1))/2;
                            head = lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation - lateral.PARENT.STATVAR.water_table_elevation;
                            contact_height = -(depth_rel2surface(1:end-1,1) - depth_rel2surface(2:end,1));
                            seepage_flux =   contact_length .* contact_height .* lateral.PARENT.ENSEMBLE{j,1}.hydraulicConductivity .* head ./ distance; %inflow!
                            
                            seepage_flux_energy = sum(seepage_flux .* lateral.PARENT.CONST.c_w .* lateral.PARENT.ENSEMBLE{j,1}.T_water,1);
                            seepage_flux = sum(seepage_flux, 1);
                            %put water in uppmost saturated cell
                            if lateral.PARENT.STATVAR.water_available
                                top_cell_saturated = find(lateral.PARENT.STATVAR.water_status(:,1)>0, 1);
                                flux(top_cell_saturated,1) = flux(top_cell_saturated,1) + seepage_flux;
                                flux_energy(top_cell_saturated,1) = flux_energy(top_cell_saturated,1) + seepage_flux_energy;
                            elseif ~isempty(lateral.PARENT.STATVAR.water_status)  %lateral.PARENT.STATVAR.water_table_top_cell > 0
                                flux(end,1) = flux(end,1) + seepage_flux;
                                flux_energy(end,1) = flux_energy(end,1) + seepage_flux_energy;
                            else
                                lateral.STATVAR.surface_runoff = lateral.STATVAR.surface_runoff + flux;
                            end
                            
                        elseif lateral.PARENT.STATVAR.water_table_elevation > lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation %outflow
                            depth_rel2surface = lateral.PARENT.STATVAR.water_table_elevation - lateral.PARENT.STATVAR.depths;
                            depth2other_worker_water_table = lateral.PARENT.STATVAR.water_table_elevation - lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation;
                            depth_rel2surface = max(0,depth_rel2surface);
                            depth_rel2surface = min(depth2other_worker_water_table,depth_rel2surface);
                            %head = (depth_rel2surface(1:end-1,1) + depth_rel2surface(2:end,1))/2;
                            head = lateral.PARENT.STATVAR.water_table_elevation - lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation;
                            contact_height = -(depth_rel2surface(1:end-1,1) - depth_rel2surface(2:end,1));
                            seepage_flux =  - contact_length .* contact_height .* lateral.PARENT.STATVAR.hydraulicConductivity .* head ./ distance; %outflow!
                            flux = flux + seepage_flux;
                            flux_energy = flux_energy + seepage_flux .* lateral.PARENT.CONST.c_w .* lateral.PARENT.STATVAR.T_water;
                        end
                    end
                end
            end
            
            
            
            %a->inflow of seepage face water
            %1. find the uppermost cell of the own worker that is saturated
            %water (water_status>0)
            %2. go down the cells of the other worker with water_status >0,
            %calculate seepage head and the flux
            %until either the bottom is reached or the cell's lower level
            %is below the water level of the own worker
            %4. add the fluxes up and add the flux to cell found in 1.
            
            %b -> outflow of seepgae face water
            %steps 2 and 3, but subrtact the water from each cell

            lateral.PARENT.STATVAR.water_flux = flux .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
            lateral.PARENT.STATVAR.water_flux_energy = flux_energy .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
            lateral.STATVAR.subsurface_run_off = sum(lateral.PARENT.STATVAR.water_flux);
            
             %modified Sep2020, moved to push
%             if ~isempty(lateral.PARENT.STATVAR.water_flux)
%                 if lateral.PARENT.STATVAR.water_table_top_cell>0
%                     lateral.PARENT.STATVAR.water_flux(lateral.PARENT.STATVAR.water_table_top_cell,1) = lateral.PARENT.STATVAR.water_flux(lateral.PARENT.STATVAR.water_table_top_cell,1) + lateral.PARENT.STATVAR.water_flux(lateral.PARENT.STATVAR.water_table_top_cell+1,1);
%                     lateral.PARENT.STATVAR.water_flux(lateral.PARENT.STATVAR.water_table_top_cell+1,:) = [];
%                     lateral.PARENT.STATVAR.water_flux_energy(lateral.PARENT.STATVAR.water_table_top_cell,1) = lateral.PARENT.STATVAR.water_flux_energy(lateral.PARENT.STATVAR.water_table_top_cell,1) + lateral.PARENT.STATVAR.water_flux_energy(lateral.PARENT.STATVAR.water_table_top_cell+1,1);
%                     lateral.PARENT.STATVAR.water_flux_energy(lateral.PARENT.STATVAR.water_table_top_cell+1,:) = [];
%                     
%                 end
%             end
            
        end

        function lateral = push(lateral, tile)
            if ~isempty(lateral.PARENT.STATVAR.water_flux)
                
                %modified Sep2020
                if lateral.PARENT.STATVAR.water_table_top_cell>0
                    lateral.PARENT.STATVAR.water_flux(lateral.PARENT.STATVAR.water_table_top_cell,1) = lateral.PARENT.STATVAR.water_flux(lateral.PARENT.STATVAR.water_table_top_cell,1) + lateral.PARENT.STATVAR.water_flux(lateral.PARENT.STATVAR.water_table_top_cell+1,1);
                    lateral.PARENT.STATVAR.water_flux(lateral.PARENT.STATVAR.water_table_top_cell+1,:) = [];
                    lateral.PARENT.STATVAR.water_flux_energy(lateral.PARENT.STATVAR.water_table_top_cell,1) = lateral.PARENT.STATVAR.water_flux_energy(lateral.PARENT.STATVAR.water_table_top_cell,1) + lateral.PARENT.STATVAR.water_flux_energy(lateral.PARENT.STATVAR.water_table_top_cell+1,1);
                    lateral.PARENT.STATVAR.water_flux_energy(lateral.PARENT.STATVAR.water_table_top_cell+1,:) = [];
                    
                end
                
                lateral.PARENT.STATVAR.water_up = 0;
                lateral.PARENT.STATVAR.water_up_energy = 0;
                
                CURRENT = lateral.PARENT.TOP.NEXT; %find correct stratigraphy class
                size_to_here = size(CURRENT.STATVAR.energy,1);
                while ~(strcmp(class(CURRENT), 'Bottom')) && size_to_here < size(lateral.PARENT.STATVAR.water_flux,1)
                    CURRENT = CURRENT.NEXT;
                    size_to_here = size_to_here + size(CURRENT.STATVAR.energy,1);
                end

                while ~(strcmp(class(CURRENT), 'Top'))
                    CURRENT = lateral3D_push_water_unconfined_aquifer(CURRENT, lateral);
                    CURRENT = CURRENT.PREVIOUS;
                end
            end
        end
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 0;
            if t + lateral.PARENT.IA_TIME_INCREMENT >= lateral.PARA.ia_time_next - 1e-7
                lateral.PARENT.ACTIVE(i,1) = 1;
                lateral.PARA.ia_time_next = t + lateral.PARENT.IA_TIME_INCREMENT + lateral.PARA.ia_time_increment;
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
            
            ground.PARA.default_value.hardBottom_cutoff = {0.03};
            ground.PARA.comment.hardBottom_cutoff = {'hard bottom  = no water flow if saturated and water content below [vol. water content, -]'};
            
            ground.PARA.default_value.ia_time_increment = {0.25};
            ground.PARA.comment.ia_time_increment ={'time step [days], must be multiple of LATERAL class timestep'};
        end
        
    end
    
end
