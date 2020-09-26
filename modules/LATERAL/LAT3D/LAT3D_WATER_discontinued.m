
classdef LAT3D_WATER < BASE_LATERAL

    
    methods
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = 24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.hardBottom_cutoff = 0.03; %hard bottom if saturated and water content below
            
            

            lateral.PARA.ia_time_increment = 0.25; %must be a multiple of the time increment of the main lateral class
            lateral.PARA.ia_time_next = [];
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.subsurface_run_off = [];
            lateral.STATVAR.surface_runoff = 0;
        end
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function lateral = finalize_init(lateral)
            lateral.STATVAR.subsurface_run_off = 0;
            lateral.STATVAR.surface_runoff = 0;
        end
        
        function lateral = pull(lateral)
            
            lateral.PARENT.STATVAR.depths = [];
            lateral.PARENT.STATVAR.water_status = [];
            lateral.PARENT.STATVAR.hydraulicConductivity = [];
            lateral.PARENT.STATVAR.water_table_elevation = [];
            lateral.PARENT.STATVAR.water_available = 0;
            lateral.PARENT.STATVAR.T_water = [];
            lateral.PARENT.STATVAR.water_table_top_cell = 0;         
            lateral.PARENT.STATVAR.aquifer_index = []; %same size as lateral.PARENT.depths - uppermost aquifer gets index 1, then 2, etc. - 0 for hard layers          
            lateral.PARENT.STATVAR.aquifer_cell_info=[];
            
            %variables required confined aquifers

            lateral.PARENT.STATVAR.head = []; %same size as number of aquifers - 1D head 
            lateral.PARENT.STATVAR.head_unknown = [];
            lateral.PARENT.STATVAR.status_head = [];
            
            lateral.PARENT.STATVAR.effective_hydraulic_conductivity =[]; %same size as number of aquifers sum(K_eff/d_eff .* A)
            lateral.PARENT.STATVAR.empty_volume_left = []; 
            lateral.PARENT.STATVAR.available_water_volume = [];
            lateral.STATVAR.global_info=cell(lateral.PARENT.PARA.num_realizations,1);
            %lateral.STATVAR.saturated = [];
            
            lateral.TEMP.open_system = 1; %start with open system
            lateral.TEMP.aquifer_index_count = 1;
            lateral.TEMP.head_space_available = 1;
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            while ~(strcmp(class(CURRENT), 'Bottom')) && lateral.TEMP.open_system == 1
                CURRENT = lateral3D_pull_water_general_aquifer(CURRENT, lateral);
                CURRENT = CURRENT.NEXT;
            end
            
            %compute empty volume and available water for each aquifer
            number_of_own_aquifers = max(lateral.PARENT.STATVAR.aquifer_index);
            for j=1:number_of_own_aquifers
                a=find(lateral.PARENT.STATVAR.aquifer_index==j, 1, 'first');
                b=find(lateral.PARENT.STATVAR.aquifer_index==j, 1, 'last');
                
                lateral.PARENT.STATVAR.empty_volume_left = [lateral.PARENT.STATVAR.empty_volume_left; -sum(lateral.PARENT.STATVAR.water_status(a:b,1).* double(lateral.PARENT.STATVAR.water_status(a:b,1)<0))];
                lateral.PARENT.STATVAR.available_water_volume = [lateral.PARENT.STATVAR.available_water_volume; sum(lateral.PARENT.STATVAR.water_status(a:b,1).* double(lateral.PARENT.STATVAR.water_status(a:b,1)>0))];
                lateral.PARENT.STATVAR.status_head = [lateral.PARENT.STATVAR.status_head; mean(lateral.PARENT.STATVAR.head_unknown(a:b,1))];
            end
            
            %make a variable in STATVAR2ALL if all aqifers are unconfined
            
        end
        
        function lateral = get_derivatives(lateral) %no need to loop through stratigraphy, all the information is in lateral.PARENT
            
            
            %calculate the bulk indices for the overlap between each pair of aquifers
            number_of_own_aquifers = max(lateral.PARENT.STATVAR.aquifer_index);
            for j=1:number_of_own_aquifers
                a=find(lateral.PARENT.STATVAR.aquifer_index==j, 1, 'first');
                b=find(lateral.PARENT.STATVAR.aquifer_index==j, 1, 'last');

                lateral.PARENT.STATVAR.aquifer_cell_info = [lateral.PARENT.STATVAR.aquifer_cell_info; [j a b]];
                cell_1 = lateral.PARENT.STATVAR.depths(a:b+1,1);
                
                for i=1:size(lateral.PARENT.ENSEMBLE,1)
                    if lateral.PARENT.PARA.connected(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{i,1}.index)
                        number_of_neighboring_aquifers = max(lateral.PARENT.ENSEMBLE{i,1}.aquifer_index);
                        
                        for k=1:number_of_neighboring_aquifers
                            c=find(lateral.PARENT.ENSEMBLE{i,1}.aquifer_index==k, 1, 'first');
                            d=find(lateral.PARENT.ENSEMBLE{i,1}.aquifer_index==k, 1, 'last');
                            cell_2 = lateral.PARENT.ENSEMBLE{i,1}.depths(c:d+1,1);
                            
                            overlap = get_overlap_aquifers(lateral.PARENT, cell_1, cell_2);
                            overlap=[overlap zeros(size(overlap,1),3)];
                            
                            
                            effective_hydraulic_conductivity=0;
                            effective_hydraulic_conductivity_times_head = 0;
                            effective_hydraulic_conductivity_times_head_ensemble = 0;
                            for l=1:size(overlap,1)
                                ehc = lateral.PARENT.STATVAR.hydraulicConductivity(overlap(l,1)+a-1,1) .* lateral.PARENT.ENSEMBLE{i,1}.hydraulicConductivity(overlap(l,2)+c-1,1)./ ...
                                    (lateral.PARENT.STATVAR.hydraulicConductivity(overlap(l,1)+a-1,1) .* lateral.PARENT.PARA.distance(lateral.PARENT.ENSEMBLE{i,1}.index, lateral.PARENT.STATVAR.index) +...
                                    lateral.PARENT.ENSEMBLE{i,1}.hydraulicConductivity(overlap(l,2)+c-1,1) .* lateral.PARENT.PARA.distance(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{i,1}.index)) .* overlap(l,3);
                                ehc(isnan(ehc)) = 0;
                                ehc = ehc .* lateral.PARENT.PARA.contact_length(lateral.PARENT.ENSEMBLE{i,1}.index, lateral.PARENT.STATVAR.index);
                                overlap(l,5) = ehc;
                                overlap(l,6) = ehc .* (double(lateral.PARENT.STATVAR.water_status(overlap(l,1)+a-1,1)>0) .* lateral.PARENT.STATVAR.head(overlap(l,1)+a-1,1) + ...
                                    double(lateral.PARENT.STATVAR.water_status(overlap(l,1)+a-1,1)<0) .* overlap(l,4));
                                overlap(l,7) = ehc .* (double(lateral.PARENT.ENSEMBLE{i,1}.water_status(overlap(l,2)+c-1,1)>0) .* lateral.PARENT.ENSEMBLE{i,1}.head(overlap(l,2)+c-1,1) + ...
                                    double(lateral.PARENT.ENSEMBLE{i,1}.water_status(overlap(l,2)+c-1,1)<0) .* overlap(l,4));
                                                                
                                effective_hydraulic_conductivity =  effective_hydraulic_conductivity + ehc;
                                effective_hydraulic_conductivity_times_head = effective_hydraulic_conductivity_times_head + overlap(l,6);
                                effective_hydraulic_conductivity_times_head_ensemble = effective_hydraulic_conductivity_times_head_ensemble + overlap(l,7);
                            end
                            %effective_hydraulic_conductivity = effective_hydraulic_conductivity.*lateral.PARENT.PARA.contact_length(lateral.PARENT.ENSEMBLE{i,1}.index, lateral.PARENT.STATVAR.index);
                            
                            lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k} = overlap; %store for later when the heads are known
                            lateral.PARENT.STATVAR.effective_hydraulic_conductivity =[lateral.PARENT.STATVAR.effective_hydraulic_conductivity; ...
                                [lateral.PARENT.ENSEMBLE{i,1}.index j k  effective_hydraulic_conductivity effective_hydraulic_conductivity_times_head effective_hydraulic_conductivity_times_head_ensemble]];
                        end
                    end
                    
                end
            end
            
            %make a column vector
            lateral.PARENT.STATVAR.effective_hydraulic_conductivity = lateral.PARENT.STATVAR.effective_hydraulic_conductivity';
            lateral.PARENT.STATVAR.effective_hydraulic_conductivity = lateral.PARENT.STATVAR.effective_hydraulic_conductivity(:);
            
            %send everything to central worker

            triangular_coupling = sum(logical(lateral.PARENT.PARA.connected(:,lateral.PARENT.PARA.central_worker)) & logical(lateral.PARENT.PARA.connected(:,lateral.PARENT.STATVAR.index)))>0 && ...
                lateral.PARENT.STATVAR.index~=lateral.PARENT.PARA.central_worker;

            if triangular_coupling || ~(lateral.PARENT.STATVAR.index == lateral.PARENT.PARA.central_worker || lateral.PARENT.PARA.connected(lateral.PARENT.STATVAR.index, lateral.PARENT.PARA.central_worker)) %own worker not connected to central worker
                %pack and send
                data_package_out = pack(lateral, {'effective_hydraulic_conductivity'; 'empty_volume_left'; 'available_water_volume'; 'status_head'});
                labSend(data_package_out, lateral.PARENT.PARA.central_worker, 1);
                
            end
            %central worker prepares own information and receives
            if lateral.PARENT.STATVAR.index == lateral.PARENT.PARA.central_worker 
                %own worker
                number_of_aquifers_per_worker = zeros(lateral.PARENT.PARA.num_realizations,1);
                %lateral.STATVAR.global_info{lateral.PARENT.STATVAR.index,1}.STATVAR.head = lateral.PARENT.STATVAR.head;
                number_of_aquifers_per_worker(lateral.PARENT.STATVAR.index,1) = size(lateral.PARENT.STATVAR.empty_volume_left,1);
                
                lateral.STATVAR.global_info{lateral.PARENT.STATVAR.index,1}.STATVAR.effective_hydraulic_conductivity = lateral.PARENT.STATVAR.effective_hydraulic_conductivity;
                lateral.STATVAR.global_info{lateral.PARENT.STATVAR.index,1}.STATVAR.empty_volume_left = lateral.PARENT.STATVAR.empty_volume_left;
                lateral.STATVAR.global_info{lateral.PARENT.STATVAR.index,1}.STATVAR.available_water_volume = lateral.PARENT.STATVAR.available_water_volume;
                lateral.STATVAR.global_info{lateral.PARENT.STATVAR.index,1}.STATVAR.status_head = lateral.PARENT.STATVAR.status_head;
                %connected workers
                for i=1:size(lateral.PARENT.ENSEMBLE,1)
                    %lateral.STATVAR.global_info{lateral.PARENT.ENSEMBLE{i,1}.index,1}.STATVAR.head = lateral.PARENT.ENSEMBLE{i,1}.head;
                    number_of_aquifers_per_worker(lateral.PARENT.ENSEMBLE{i,1}.index,1) = size(lateral.PARENT.ENSEMBLE{i,1}.empty_volume_left,1);
                    lateral.STATVAR.global_info{lateral.PARENT.ENSEMBLE{i,1}.index,1}.STATVAR.effective_hydraulic_conductivity = []; %info already provided by central_worker
                    lateral.STATVAR.global_info{lateral.PARENT.ENSEMBLE{i,1}.index,1}.STATVAR.empty_volume_left = lateral.PARENT.ENSEMBLE{i,1}.empty_volume_left;
                    lateral.STATVAR.global_info{lateral.PARENT.ENSEMBLE{i,1}.index,1}.STATVAR.available_water_volume = lateral.PARENT.ENSEMBLE{i,1}.available_water_volume;
                    lateral.STATVAR.global_info{lateral.PARENT.ENSEMBLE{i,1}.index,1}.STATVAR.status_head = lateral.PARENT.ENSEMBLE{i,1}.status_head;
                end
                
                %only the workers 
                sending_workers = ~lateral.PARENT.PARA.connected(lateral.PARENT.STATVAR.index,:);
                not_sending_workers = find(sending_workers == 0);
                for i=1:size(not_sending_workers,2)
                    sending_workers(not_sending_workers(1,i)) = ...                    
                          sum(logical(lateral.PARENT.PARA.connected(:,not_sending_workers(1,i))) & logical(lateral.PARENT.PARA.connected(:,lateral.PARENT.STATVAR.index)))>0;
                end
                sending_workers(1,lateral.PARENT.STATVAR.index) = 0;
                sending_workers = find(sending_workers == 1);
                
                %sending_workers
                for i=sending_workers
                    %receive and unpack
                    data_package_in = labReceive(i, 1);
                    %and unpack
                    lateral = unpack(lateral, data_package_in, i);
                    number_of_aquifers_per_worker(i,1) = size(lateral.STATVAR.global_info{i,1}.STATVAR.empty_volume_left,1);
                end
                
                total_number_of_aquifers = sum(number_of_aquifers_per_worker,1);
                
                K_prime = zeros(total_number_of_aquifers, total_number_of_aquifers);
                K_prime_times_head = zeros(total_number_of_aquifers, total_number_of_aquifers);
                head_known =  zeros(total_number_of_aquifers,1);
                
                offset=0;
                offset_store=0;
                for i = 1:size(lateral.STATVAR.global_info,1)
                    %head(offset+1:offset+size(lateral.STATVAR.global_info{i,1}.STATVAR.head,1),1) = lateral.STATVAR.global_info{i,1}.STATVAR.head;
                    head_known(offset+1:offset+size(lateral.STATVAR.global_info{i,1}.STATVAR.empty_volume_left,1),1) = lateral.STATVAR.global_info{i,1}.STATVAR.status_head < 2;
                    %head_known(offset+1,1) = 1;  %chnage this, this variable must be created at class level, using the open_system variable!
                    offset = offset + size(lateral.STATVAR.global_info{i,1}.STATVAR.empty_volume_left,1);
                    offset_store=[offset_store; offset];
                end
                
                for i = 1:size(lateral.STATVAR.global_info,1)
                    for j=1:6:size(lateral.STATVAR.global_info{i,1}.STATVAR.effective_hydraulic_conductivity,1)
                        i1 = offset_store(i,1) + lateral.STATVAR.global_info{i,1}.STATVAR.effective_hydraulic_conductivity(j+2-1,1);
                        i2 = offset_store(lateral.STATVAR.global_info{i,1}.STATVAR.effective_hydraulic_conductivity(j+1-1,1),1) + lateral.STATVAR.global_info{i,1}.STATVAR.effective_hydraulic_conductivity(j+3-1,1);
                        K_prime(i1, i2) = lateral.STATVAR.global_info{i,1}.STATVAR.effective_hydraulic_conductivity(j+4-1,1);
                        K_prime(i2, i1) = K_prime(i1, i2);
                        K_prime_times_head(i1, i2) = lateral.STATVAR.global_info{i,1}.STATVAR.effective_hydraulic_conductivity(j+5-1,1);
                        K_prime_times_head(i2, i1) = lateral.STATVAR.global_info{i,1}.STATVAR.effective_hydraulic_conductivity(j+6-1,1);
                    end
                end
                %----delete later---
                lateral.STATVAR.K_prime=K_prime;
                lateral.STATVAR.K_prime_times_head=K_prime_times_head;
                lateral.STATVAR.head_known=head_known;
                recomputed_head = head_known.* 0;
                %-------------
                
                enumeration_pressure_unknown = find(~head_known);
                enumeration_pressure_known = find(head_known);
                lin_matrix = zeros(sum(double(~head_known)), sum(double(~head_known)));
                
                for i=1:size(enumeration_pressure_unknown,1)
                    for j=1:i
                        if i==j
                            lin_matrix(i,i) = sum(K_prime(:,enumeration_pressure_unknown(i,1)),1);
                        else
                            lin_matrix(i,j) = -K_prime(enumeration_pressure_unknown(i,1), enumeration_pressure_unknown(j,1));
                            lin_matrix(j,i) =lin_matrix(i,j);
                        end
                    end
                end
                
                lin_vector = zeros(sum(double(~head_known)),1);
                for i=1:size(enumeration_pressure_unknown,1)
                    for j=1:size(enumeration_pressure_known,1)
                        lin_vector(i,1) = lin_vector(i,1) + K_prime_times_head(enumeration_pressure_known(j,1), enumeration_pressure_unknown(i,1));
                    end
                end
                
                %head(~head_known) = linsolve(lin_matrix,lin_vector);
                
                unknown_heads = linsolve(lin_matrix,lin_vector);
                
                recomputed_head(~head_known) = unknown_heads;
                
                lateral.PARENT.STATVAR.recomputed_head=recomputed_head;
                
                for i=1:size(enumeration_pressure_unknown,1)
                    K_prime_times_head(enumeration_pressure_unknown(i,1),:) = K_prime(enumeration_pressure_unknown(i,1),:) .* unknown_heads(i,1);
                end
                
                %flux_matrix = K_prime.*(head-head');
                flux_matrix = K_prime_times_head - K_prime_times_head';
                
                lateral.STATVAR.flux_matrix = flux_matrix;
                
                total_fluxes = sum(flux_matrix,1)';
                
                fluxes = {};
                offset = 0;
                for i = 1:lateral.PARENT.PARA.num_realizations
                    fluxes{i,1} = total_fluxes(offset+1:offset+number_of_aquifers_per_worker(i,1),1);
                    offset = offset + number_of_aquifers_per_worker(i,1);
                end

                lateral.STATVAR.fluxes = fluxes;
%                 
                %pack and send back to all workers
                for i= 1:lateral.PARENT.PARA.num_realizations
                    if i ~=lateral.PARENT.STATVAR.index
                        data_package_out = pack(lateral, {'recomputed_head'});
                        labSend(data_package_out, i, 1);
                    end
                end

            end
            %labBarrier();

            
            %receive and unpack
            if lateral.PARENT.STATVAR.index ~= lateral.PARENT.PARA.central_worker
                data_package_in = labReceive(lateral.PARENT.PARA.central_worker, 1);
                lateral = unpack(lateral, data_package_in, i);
            end

            

            
            %process
            
            
            
            
            labBarrier();
            
            
            
            %determine the heads
            
            %send everything back
            
%             %calculate fluxes
%             flux = lateral.PARENT.STATVAR.hydraulicConductivity .* 0;
%             flux_energy = flux;
%             for j=1:size(lateral.PARENT.ENSEMBLE,1)
%                 if lateral.PARENT.PARA.connected(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index)  %add if water is available at all
%                     %flow between saturated cells
%                     contact_length = lateral.PARENT.PARA.contact_length(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
%                     distance = lateral.PARENT.PARA.distance(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
%                     for i=1:size(lateral.PARENT.ENSEMBLE{j,1}.overlap,1) %does not do anything when overlap is empty!
%                         
%                         hc1 = lateral.PARENT.STATVAR.hydraulicConductivity(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1);
%                         hc2 = lateral.PARENT.ENSEMBLE{j,1}.hydraulicConductivity(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,2),1);
%                         hydraulicConductivity = hc1 .* hc2 ./ (hc1 .* distance./2 + hc2 .* distance./2); %change to different distances laters
%                         hydraulicConductivity(isnan(hydraulicConductivity)) = 0;
% 
%                         flux_i = contact_length .* lateral.PARENT.ENSEMBLE{j,1}.overlap(i,3) .* hydraulicConductivity .* ...
%                             -(lateral.PARENT.STATVAR.water_table_elevation - lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation);
%                         flux(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) = flux(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) + flux_i;
%                         flux_energy(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) = flux_energy(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) + flux_i .* lateral.PARENT.CONST.c_w .* ...
%                             (double(flux_i<0).*lateral.PARENT.STATVAR.T_water(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) + ...
%                             double(flux_i>=0) .* lateral.PARENT.ENSEMBLE{j,1}.T_water(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,2),1));
%                         
%                     end
%                     %seepage flow of saturated cells above other cell's water table
%                     if lateral.PARENT.STATVAR.water_available || lateral.PARENT.ENSEMBLE{j,1}.water_available
%                         if lateral.PARENT.STATVAR.water_table_elevation < lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation %inflow
%                             depth_rel2surface = lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation - lateral.PARENT.ENSEMBLE{j,1}.depths;
%                             depth2other_worker_water_table = lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation - lateral.PARENT.STATVAR.water_table_elevation;
%                             depth_rel2surface = max(0,depth_rel2surface);
%                             depth_rel2surface = min(depth2other_worker_water_table,depth_rel2surface);
%                             head = (depth_rel2surface(1:end-1,1) + depth_rel2surface(2:end,1))/2;
%                             contact_height = -(depth_rel2surface(1:end-1,1) - depth_rel2surface(2:end,1));
%                             seepage_flux =   contact_length .* contact_height .* lateral.PARENT.ENSEMBLE{j,1}.hydraulicConductivity .* head ./ distance; %inflow!
%                             
%                             seepage_flux_energy = sum(seepage_flux .* lateral.PARENT.CONST.c_w .* lateral.PARENT.ENSEMBLE{j,1}.T_water,1);
%                             seepage_flux = sum(seepage_flux, 1);
%                             %put water in uppmost saturated cell
%                             if lateral.PARENT.STATVAR.water_available
%                                 top_cell_saturated = find(lateral.PARENT.STATVAR.water_status(:,1)>0, 1);
%                                 flux(top_cell_saturated,1) = flux(top_cell_saturated,1) + seepage_flux;
%                                 flux_energy(top_cell_saturated,1) = flux_energy(top_cell_saturated,1) + seepage_flux_energy;
%                             elseif ~isempty(lateral.PARENT.STATVAR.water_status)  %lateral.PARENT.STATVAR.water_table_top_cell > 0
%                                 flux(end,1) = flux(end,1) + seepage_flux;
%                                 flux_energy(end,1) = flux_energy(end,1) + seepage_flux_energy;
%                             else
%                                 lateral.STATVAR.surface_runoff = lateral.STATVAR.surface_runoff + flux;
%                             end
%                             
%                         elseif lateral.PARENT.STATVAR.water_table_elevation > lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation %outflow
%                             depth_rel2surface = lateral.PARENT.STATVAR.water_table_elevation - lateral.PARENT.STATVAR.depths;
%                             depth2other_worker_water_table = lateral.PARENT.STATVAR.water_table_elevation - lateral.PARENT.ENSEMBLE{j,1}.water_table_elevation;
%                             depth_rel2surface = max(0,depth_rel2surface);
%                             depth_rel2surface = min(depth2other_worker_water_table,depth_rel2surface);
%                             head = (depth_rel2surface(1:end-1,1) + depth_rel2surface(2:end,1))/2;
%                             contact_height = -(depth_rel2surface(1:end-1,1) - depth_rel2surface(2:end,1));
%                             seepage_flux =  - contact_length .* contact_height .* lateral.PARENT.STATVAR.hydraulicConductivity .* head ./ distance; %outflow!
%                             flux = flux + seepage_flux;
%                             flux_energy = flux_energy + seepage_flux .* lateral.PARENT.CONST.c_w .* lateral.PARENT.STATVAR.T_water;
%                         end
%                     end
%                 end
%             end
%             
%             
%             
%             %a->inflow of seepage face water
%             %1. find the uppermost cell of the own worker that is saturated
%             %water (water_status>0)
%             %2. go down the cells of the other worker with water_status >0,
%             %calculate seepage head and the flux
%             %until either the bottom is reached or the cell's lower level
%             %is below the water level of the own worker
%             %4. add the fluxes up and add the flux to cell found in 1.
%             
%             %b -> outflow of seepgae face water
%             %steps 2 and 3, but subrtact the water from each cell
% 
%             lateral.PARENT.STATVAR.water_flux = flux .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
%             lateral.PARENT.STATVAR.water_flux_energy = flux_energy .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
            
            
        end

        
        
        function lateral = push(lateral, forcing)
            if ~isempty(lateral.PARENT.STATVAR.water_flux)
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
        
        
        function data_package = pack(lateral, variable_list) %transform lateral.STATVAR into column vector ready to send
            data_package = [];
            for i=1:size(variable_list,1)
                data_package=[data_package; size(variable_list{i,1},2); double(variable_list{i,1})']; % # of characters followed by characters as doubles
                data_package=[data_package; size(lateral.PARENT.STATVAR.(variable_list{i,1}),1); lateral.PARENT.STATVAR.(variable_list{i,1})]; % # of entries followed by values
            end
        end
        
        function lateral = unpack(lateral, data_package, index) %read received column vector and transform into STATVAR
            i=1;
            while i<=size(data_package,1)
               variable_name = char(data_package(i+1:i+data_package(i,1),1)');
               i = i + data_package(i,1) + 1;
               STATVAR.(variable_name) = data_package(i+1:i+data_package(i,1),1);
               i = i + data_package(i,1) + 1;
            end
            lateral.STATVAR.global_info{index,1}.STATVAR = STATVAR;
            %number_of_aquifers = size(STATVAR.head,1);
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


