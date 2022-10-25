%========================================================================
% CryoGrid LATERAL_IA class LAT3D_WATER_UNCONFINED_AQUIFER_OVERLAND_FLOW 
% simulates lateral water flow between pairs of CryoGrid stratigraphies
% for the topmost unconfined aquifer, as well as overland flow according to
% Gauckler-Manning equation
% S. Westermann, Oct 2020, June 2021
%========================================================================


classdef LAT3D_WATER_UNCONFINED_AQUIFER_OVERLAND_FLOW < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = []; %24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.hardBottom_cutoff = []; %0.03; %hard bottom if saturated and water content below
            lateral.PARA.GaMa_coefficient = []; %Gauckler-Manning coefficient, https://en.wikipedia.org/wiki/Manning_formula
            lateral.PARA.tortuosity = 1;
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
            lateral.PARENT.STATVAR.water_depth = 0; %depth of the free water at the surface flowing according to Gauckler-Manning 
            lateral.PARENT.STATVAR.ground_surface_elevation = []; %NEW SEBASTIAN for hillslope flow between grid cells

            lateral.TEMP.open_system = 1; %start with open system

            CURRENT = lateral.PARENT.TOP.NEXT;
            CURRENT = lateral3D_pull_water_overland_flow(CURRENT, lateral);
            
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
                            head = (depth_rel2surface(1:end-1,1) + depth_rel2surface(2:end,1))/2;
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
                            head = (depth_rel2surface(1:end-1,1) + depth_rel2surface(2:end,1))/2;
                            contact_height = -(depth_rel2surface(1:end-1,1) - depth_rel2surface(2:end,1));
                            seepage_flux =  - contact_length .* contact_height .* lateral.PARENT.STATVAR.hydraulicConductivity .* head ./ distance; %outflow!
                            flux = flux + seepage_flux;
                            flux_energy = flux_energy + seepage_flux .* lateral.PARENT.CONST.c_w .* lateral.PARENT.STATVAR.T_water;
                        end
                    end
                    % overland flow according to Gauckler-Manning
                    if (lateral.PARENT.STATVAR.water_depth > 1e-4 && lateral.PARENT.ENSEMBLE{j,1}.water_depth > 1e-4  && ~isempty(lateral.PARENT.STATVAR.T_water) && ~isempty(lateral.PARENT.ENSEMBLE{j,1}.T_water)) ...
                        || (lateral.PARENT.STATVAR.water_depth <= 1e-4 && lateral.PARENT.ENSEMBLE{j,1}.water_depth > 1e-4  && ~isempty(lateral.PARENT.ENSEMBLE{j,1}.T_water)  && ~isempty(lateral.PARENT.STATVAR.T_water) && lateral.PARENT.STATVAR.depths(1,1) < lateral.PARENT.ENSEMBLE{j,1}.depths(1,1)) ...
                        || (lateral.PARENT.STATVAR.water_depth > 1e-4 && lateral.PARENT.ENSEMBLE{j,1}.water_depth <= 1e-4  && ~isempty(lateral.PARENT.ENSEMBLE{j,1}.T_water)  && ~isempty(lateral.PARENT.STATVAR.T_water) && lateral.PARENT.STATVAR.depths(1,1) > lateral.PARENT.ENSEMBLE{j,1}.depths(1,1)) 
                        if lateral.PARENT.STATVAR.depths(1,1) < lateral.PARENT.ENSEMBLE{j,1}.depths(1,1) && lateral.PARENT.STATVAR.water_depth < 1e-4
                            lateral.PARENT.STATVAR.water_depth = 1e-4;
                            lateral.PARENT.STATVAR.max_flow = 1e-2 .* lateral.PARENT.STATVAR.area_flow;
                        end
                        if lateral.PARENT.STATVAR.depths(1,1) > lateral.PARENT.ENSEMBLE{j,1}.depths(1,1) && lateral.PARENT.ENSEMBLE{j,1}.water_depth < 1e-4
                            lateral.PARENT.ENSEMBLE{j,1}.water_depth = 1e-4;
                            lateral.PARENT.ENSEMBLE{j,1}.max_flow = 1e-2 .* lateral.PARENT.ENSEMBLE{j,1}.area_flow;
                        end
                        
                        
                        %could be necessary to  limit flow for numerical
                        %stability?
                        gradient = -(lateral.PARENT.STATVAR.depths(1,1) - lateral.PARENT.ENSEMBLE{j,1}.depths(1,1)) ./ (distance.*lateral.PARA.tortuosity);
                        if gradient < 0  %own realization higher
                            T_water =  lateral.PARENT.STATVAR.T_water(1,1);
                            depth1 = min(lateral.PARENT.STATVAR.depths(1,1) - (lateral.PARENT.ENSEMBLE{j,1}.depths(1,1) - lateral.PARENT.ENSEMBLE{j,1}.water_depth), lateral.PARENT.STATVAR.water_depth);
                            depth2 = lateral.PARENT.ENSEMBLE{j,1}.water_depth;
                        elseif gradient > 0  %ensemble realization higher
                            T_water = lateral.PARENT.ENSEMBLE{j,1}.T_water(1,1);
                            depth1 = min(lateral.PARENT.ENSEMBLE{j,1}.depths(1,1) - (lateral.PARENT.STATVAR.depths(1,1) - lateral.PARENT.STATVAR.water_depth), lateral.PARENT.ENSEMBLE{j,1}.water_depth);
                            depth2 = lateral.PARENT.STATVAR.water_depth;
                        end
                        
                        
                        %velocity = lateral.PARA.GaMa_coefficient .* real(min(lateral.PARENT.STATVAR.water_depth.^(2/3), lateral.PARENT.ENSEMBLE{j,1}.water_depth.^(2/3)) .* abs(gradient).^0.5);
                        velocity = lateral.PARA.GaMa_coefficient .* real(((depth1+depth2)./2).^(2/3) .* abs(gradient).^0.5);
                        %flow =  velocity .* min(lateral.PARENT.STATVAR.water_depth, lateral.PARENT.ENSEMBLE{j,1}.water_depth) .* contact_length; 
                        flow =  velocity .* (depth1+depth2)./2 .* contact_length; 
%                         if flow <0
%                             flow = max(flow, -lateral.PARENT.STATVAR.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec));
%                         end
                        %flow that would lead to the same water level in
                        %both tiles
                        max_flow2same_level = lateral.PARENT.STATVAR.area_flow .* lateral.PARENT.ENSEMBLE{j,1}.area_flow .* (lateral.PARENT.STATVAR.depths(1,1) - lateral.PARENT.ENSEMBLE{j,1}.depths(1,1)) ./ (lateral.PARENT.STATVAR.area_flow + lateral.PARENT.ENSEMBLE{j,1}.area_flow) ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec);
                        
                        %TEST Sebastian                        
                        %flow = sign(gradient) .* min(abs(max_flow2same_level./8), min(min(flow, lateral.PARENT.STATVAR.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec)), lateral.PARENT.ENSEMBLE{j,1}.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec)));


                        if gradient < 0 %own realization higher, looses water
                            flow = sign(gradient) .* min(abs(max_flow2same_level./8), min(flow, lateral.PARENT.STATVAR.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec)));
                        else
                            flow = sign(gradient) .* min(abs(max_flow2same_level./8), min(flow, lateral.PARENT.ENSEMBLE{j,1}.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec)));
                        end
                    %end TEST Sebastian
                    
                    
                        flow_energy = flow .* lateral.PARENT.CONST.c_w .* T_water;
                        
                        flux(1,1) = flux(1,1) + flow;
                        flux_energy(1,1) = flux_energy(1,1) + flow_energy;
                        
                        lateral.STATVAR.flow = flow;
                        lateral.STATVAR.max_flow = -lateral.PARENT.STATVAR.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec);
                        lateral.STATVAR.water_depth = lateral.PARENT.STATVAR.water_depth;
                        
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
            
            ground.PARA.default_value.hardBottom_cutoff = {0.03};
            ground.PARA.comment.hardBottom_cutoff = {'hard bottom  = no water flow if saturated and water content below [vol. water content, -]'};
            
            ground.PARA.default_value.GaMa_coefficient = {15};
            ground.PARA.comment.GaMa_coefficient = {'Gauckler-Manning coefficient, https://en.wikipedia.org/wiki/Manning_formula'};
            
            ground.PARA.default_value.tortuosity = {1};
            ground.PARA.comment.tortuosity = {'multiply direct distance with this factor to get flow path length'};
            
            ground.PARA.default_value.ia_time_increment = {0.25};
            ground.PARA.comment.ia_time_increment ={'time step [days], must be multiple of LATERAL class timestep'};
        end
    end
    
end