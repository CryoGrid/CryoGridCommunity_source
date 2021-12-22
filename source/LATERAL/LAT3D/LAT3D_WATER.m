%========================================================================
% CryoGrid LATERAL_IA class LAT3D_WATER
% simulates lateral water flow between pairs of CryoGrid stratigraphies.
% NOTE: does not work stable yet, DO NOT USE!
% S. Westermann, Oct 2020
%========================================================================


classdef LAT3D_WATER < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = []; %24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.hardBottom_cutoff = []; %0.03; %hard bottom if saturated and water content below
            lateral.PARA.ia_time_increment = []; %0.25; %must be a multiple of the time increment of the main lateral class
            %lateral.PARA.ia_time_next = [];
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.subsurface_run_off = [];
            lateral.STATVAR.surface_runoff = []; %0;
        end

        function lateral = finalize_init(lateral, tile)
            lateral.STATVAR.subsurface_run_off = 0;
            lateral.STATVAR.surface_runoff = 0;
        end
        
        %---- time integration-----------------
        
        function lateral = pull(lateral, tile)
            
            lateral.PARENT.STATVAR.depths = [];
            lateral.PARENT.STATVAR.water_status = [];
            lateral.PARENT.STATVAR.hydraulicConductivity = [];
            lateral.PARENT.STATVAR.water_table_elevation = [];
            lateral.PARENT.STATVAR.water_available = 0;
            lateral.PARENT.STATVAR.T_water = [];
            lateral.PARENT.STATVAR.water_table_top_cell = 0;         
            lateral.PARENT.STATVAR.aquifer_index = []; %same size as lateral.PARENT.depths - uppermost aquifer gets index 1, then 2, etc. - 0 for hard layers          
            lateral.STATVAR.aquifer_cell_info=[];
            lateral.PARENT.STATVAR.external_fluxes = []; %for 1D lateral classes like seepage face and reservoir
            
            %variables required confined aquifers

            lateral.PARENT.STATVAR.head = []; %same size as number of aquifers - 1D head 
            lateral.PARENT.STATVAR.head_unknown = [];
            lateral.PARENT.STATVAR.status_head = [];
            
            lateral.PARENT.STATVAR.effective_hydraulic_conductivity =[]; %same size as number of aquifers sum(K_eff/d_eff .* A)
            lateral.PARENT.STATVAR.empty_volume_left = []; 
            lateral.PARENT.STATVAR.available_water_volume = [];
            lateral.STATVAR.global_info=cell(lateral.PARENT.PARA.num_realizations,1);
            
            lateral.PARENT.STATVAR_PRIVATE.overlap_info = {};
            lateral.PARENT.STATVAR_PRIVATE.overlap_aquifer_index = [];
            
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
        
        function lateral = get_derivatives(lateral, tile) %no need to loop through stratigraphy, all the information is in lateral.PARENT
            
            %calculate the bulk indices for the overlap between each pair of aquifers
            number_of_own_aquifers = max(lateral.PARENT.STATVAR.aquifer_index);
            for j=1:number_of_own_aquifers
                a=find(lateral.PARENT.STATVAR.aquifer_index==j, 1, 'first');
                b=find(lateral.PARENT.STATVAR.aquifer_index==j, 1, 'last');

                lateral.STATVAR.aquifer_cell_info = [lateral.STATVAR.aquifer_cell_info; [j a b]];
                cell_1 = lateral.PARENT.STATVAR.depths(a:b+1,1);
                
                for i=1:size(lateral.PARENT.ENSEMBLE,1)
                    if lateral.PARENT.PARA.connected(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{i,1}.index)
                        number_of_neighboring_aquifers = max(lateral.PARENT.ENSEMBLE{i,1}.aquifer_index);
                        lateral.PARENT.ENSEMBLE{i,1}.aquifer_cell_info = [];
                        
                        for k=1:number_of_neighboring_aquifers
                            c=find(lateral.PARENT.ENSEMBLE{i,1}.aquifer_index==k, 1, 'first');
                            d=find(lateral.PARENT.ENSEMBLE{i,1}.aquifer_index==k, 1, 'last');
                            lateral.PARENT.ENSEMBLE{i,1}.aquifer_cell_info = [lateral.PARENT.ENSEMBLE{i,1}.aquifer_cell_info; [k c d]];
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
                                overlap(l,6) = ehc .* (double(lateral.PARENT.STATVAR.water_status(overlap(l,1)+a-1,1)>=0) .* lateral.PARENT.STATVAR.head(overlap(l,1)+a-1,1) + ...
                                    double(lateral.PARENT.STATVAR.water_status(overlap(l,1)+a-1,1)<0) .* overlap(l,4));
                                overlap(l,7) = ehc .* (double(lateral.PARENT.ENSEMBLE{i,1}.water_status(overlap(l,2)+c-1,1)>=0) .* lateral.PARENT.ENSEMBLE{i,1}.head(overlap(l,2)+c-1,1) + ...
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
            
            %synchronize information between workers
            %send
            data_package_out = pack(lateral, {'effective_hydraulic_conductivity'; 'empty_volume_left'; 'available_water_volume'; 'status_head'; 'external_fluxes'});
            for i=1:lateral.PARENT.PARA.num_realizations
                if i~=lateral.PARENT.STATVAR.index          
                    labSend(data_package_out, i, 1);
                end
            end
            
            labBarrier();
            
            number_of_aquifers_per_worker = zeros(lateral.PARENT.PARA.num_realizations,1);
            number_of_aquifers_per_worker(lateral.PARENT.STATVAR.index,1) = size(lateral.PARENT.STATVAR.empty_volume_left,1);
            
            lateral.STATVAR.global_info{lateral.PARENT.STATVAR.index,1}.STATVAR.effective_hydraulic_conductivity = lateral.PARENT.STATVAR.effective_hydraulic_conductivity;
            lateral.STATVAR.global_info{lateral.PARENT.STATVAR.index,1}.STATVAR.empty_volume_left = lateral.PARENT.STATVAR.empty_volume_left;
            lateral.STATVAR.global_info{lateral.PARENT.STATVAR.index,1}.STATVAR.available_water_volume = lateral.PARENT.STATVAR.available_water_volume;
            lateral.STATVAR.global_info{lateral.PARENT.STATVAR.index,1}.STATVAR.status_head = lateral.PARENT.STATVAR.status_head;
            lateral.STATVAR.global_info{lateral.PARENT.STATVAR.index,1}.STATVAR.external_fluxes = lateral.PARENT.STATVAR.external_fluxes;
            
            for i=1:lateral.PARENT.PARA.num_realizations
                if i~=lateral.PARENT.STATVAR.index
                    data_package_in = labReceive(i, 1);
                    lateral = unpack(lateral, data_package_in, i);
                    number_of_aquifers_per_worker(i,1) = size(lateral.STATVAR.global_info{i,1}.STATVAR.empty_volume_left,1);
                end
            end
            
            %assign matrices for transitions between each aquifer
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
            
            
            %solve for the aquifers with unknown heads - check if the
            %computed head is lower than the original head - in this case,
            %water will drain and create air space - set original head as
            %known in this case
            
            recompute = 1;
            while recompute
                
                enumeration_pressure_unknown = find(~head_known);
                
                %add fluxes from/to external reservoirs/seepage faces
                external_head_vector = zeros(total_number_of_aquifers,1);
                external_conductivity_vector = zeros(total_number_of_aquifers,1);
                offset = 0;
                for i=1:size(lateral.STATVAR.global_info,1)
                    for j=1:4:size(lateral.STATVAR.global_info{i,1}.STATVAR.external_fluxes,1)
                        external_head_vector(offset + lateral.STATVAR.global_info{i,1}.STATVAR.external_fluxes(j-1+1,1),1) = ...
                            external_head_vector(offset + lateral.STATVAR.global_info{i,1}.STATVAR.external_fluxes(j-1+1,1),1) + lateral.STATVAR.global_info{i,1}.STATVAR.external_fluxes(j-1+4,1);
                        external_conductivity_vector(offset + lateral.STATVAR.global_info{i,1}.STATVAR.external_fluxes(j-1+1,1),1) = ...
                            external_conductivity_vector(offset + lateral.STATVAR.global_info{i,1}.STATVAR.external_fluxes(j-1+1,1),1) + lateral.STATVAR.global_info{i,1}.STATVAR.external_fluxes(j-1+2,1);
                    end
                    offset = offset + number_of_aquifers_per_worker(i,1);
                end
                %external_conductivity_vector
                %external_head_vector
                
                %remove unconnected aquifers
                for i=1:size(enumeration_pressure_unknown,1)
                    if (sum(K_prime(:,enumeration_pressure_unknown(i,1)),1) + external_conductivity_vector(enumeration_pressure_unknown(i,1),1))==0
                        head_known(enumeration_pressure_unknown(i,1),1) = 1;
                    end
                end
                
                enumeration_pressure_unknown = find(~head_known);
                enumeration_pressure_known = find(head_known);
                lin_matrix = zeros(sum(double(~head_known)), sum(double(~head_known)));
                
                for i=1:size(enumeration_pressure_unknown,1)
                    for j=1:i
                        if i==j
                            lin_matrix(i,i) = sum(K_prime(:,enumeration_pressure_unknown(i,1)),1) + external_conductivity_vector(enumeration_pressure_unknown(i,1),1);
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
                    
                    lin_vector(i,1) = lin_vector(i,1) + external_head_vector(enumeration_pressure_unknown(i,1),1);  %add external couplings for all aquifers with unknown heads
                end
                
                %solving laplace equation with known heads as BC
                unknown_heads = linsolve(lin_matrix,lin_vector);
                if sum(isnan(unknown_heads))>0 || sum(unknown_heads==Inf)>0 || sum(unknown_heads==-Inf)>0
                        'Hallo'
                      %  unknown_heads
%                     external_head_vector(1:4)
%                     external_conductivity_vector(1:4)
                end
                
                recompute = 0;
                %check if original head is higher than recomputed one
                for i=1:size(enumeration_pressure_unknown,1)

                    if isnan(unknown_heads(i,1)) || unknown_heads(i,1)==Inf || unknown_heads(i,1)==-Inf || sum(K_prime_times_head(enumeration_pressure_unknown(i,1),:),2) >  sum(K_prime(enumeration_pressure_unknown(i,1),:) .* unknown_heads(i,1),2) 
                       head_known(enumeration_pressure_unknown(i,1),1) = 1;
                       recompute = 1;
                    end
                end

            end
            
            recomputed_head(~head_known) = unknown_heads;
            lateral.PARENT.STATVAR.recomputed_head=recomputed_head;
            
            
            %------
            %compute fluxes on cell by cell basis
            
            %flux = lateral.PARENT.STATVAR.aquifer_index.*0;
            flux_in = lateral.PARENT.STATVAR.aquifer_index.*0;
            flux_out = flux_in;
            flux_energy_in = flux_in;
            flux_energy_out = flux_in;
            
            for j=1:number_of_own_aquifers  %same as size(lateral.PARENT.ENSEMBLE{i,1}.overlap_info,1)
                global_index_own_aquifer = sum(number_of_aquifers_per_worker(1:lateral.PARENT.STATVAR.index-1,1),1) + j;
                 %lateral.STATVAR.aquifer_cell_info = [lateral.STATVAR.aquifer_cell_info; [j a b]];
                 
                for i=1:size(lateral.PARENT.ENSEMBLE,1)

                    for k=1:size(lateral.PARENT.ENSEMBLE{i,1}.overlap_info,2)

                        global_index_ensemble_aquifer =  sum(number_of_aquifers_per_worker(1:lateral.PARENT.ENSEMBLE{i,1}.index - 1,1),1) + k;
                                
                        if ~isempty(lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k}) && ~head_known(global_index_own_aquifer,1)
                            lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k}(:,6) = lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k}(:,5) .* recomputed_head(global_index_own_aquifer,1);
                        end
                        if ~isempty(lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k}) && ~head_known(global_index_ensemble_aquifer,1)
                            lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k}(:,7) = lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k}(:,5) .* recomputed_head(global_index_ensemble_aquifer,1);
                        end
                        if ~isempty(lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k})
                            
                            flux_between_aquifers = lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k}(:,7) - lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k}(:,6);
                            
                            range_own_aquifer = lateral.STATVAR.aquifer_cell_info(j,2)-1 + lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k}(:,1);
                            range_ensemble_aquifer = lateral.PARENT.ENSEMBLE{i,1}.aquifer_cell_info(k,2)-1 + lateral.PARENT.ENSEMBLE{i,1}.overlap_info{j,k}(:,2);
                            for n=1:size(flux_between_aquifers,1)
                                
                                flux_in(range_own_aquifer(n,1),1) = flux_in(range_own_aquifer(n,1),1) + double(flux_between_aquifers(n,1) > 0) .* flux_between_aquifers(n,1);
                                flux_out(range_own_aquifer(n,1),1) = flux_out(range_own_aquifer(n,1),1) - double(flux_between_aquifers(n,1) < 0) .* flux_between_aquifers(n,1);
                                
                                T_water = double(flux_between_aquifers(n,1) > 0) .* lateral.PARENT.ENSEMBLE{i,1}.T_water(range_ensemble_aquifer(n,1),1) + ...
                                    double(flux_between_aquifers(n,1) < 0) .* lateral.PARENT.STATVAR.T_water(range_own_aquifer(n,1),1);
                                flux_energy_in(range_own_aquifer(n,1),1) = flux_energy_in(range_own_aquifer(n,1),1) + double(flux_between_aquifers(n,1) > 0) .* flux_between_aquifers(n,1) .* ...
                                    (double(T_water>0) .* lateral.PARENT.CONST.c_w + double(T_water<0) .* lateral.PARENT.CONST.c_i) .* T_water;
                                flux_energy_out(range_own_aquifer(n,1),1) = flux_energy_out(range_own_aquifer(n,1),1) - double(flux_between_aquifers(n,1) < 0) .* flux_between_aquifers(n,1) .* ...
                                    (double(T_water>0) .* lateral.PARENT.CONST.c_w + double(T_water<0) .* lateral.PARENT.CONST.c_i) .* T_water;
                            end
                        end
                    end
                end
            end
            
            %external reservoirs
            for i=1:size(lateral.PARENT.STATVAR_PRIVATE.overlap_aquifer_index,1)
                aquifer_index = lateral.PARENT.STATVAR_PRIVATE.overlap_aquifer_index(i,1);
                global_index_own_aquifer = sum(number_of_aquifers_per_worker(1:lateral.PARENT.STATVAR.index-1,1),1) + aquifer_index;
                
                if ~isempty(lateral.PARENT.STATVAR_PRIVATE.overlap_info{i,1}) && ~head_known(global_index_own_aquifer,1)
                    lateral.PARENT.STATVAR_PRIVATE.overlap_info{i,1}(:,6) = lateral.PARENT.STATVAR_PRIVATE.overlap_info{i,1}(:,5) .* recomputed_head(global_index_own_aquifer,1);
                end
                
                flux_between_aquifers = lateral.PARENT.STATVAR_PRIVATE.overlap_info{i,1}(:,7) - lateral.PARENT.STATVAR_PRIVATE.overlap_info{i,1}(:,6);
                range_own_aquifer = lateral.STATVAR.aquifer_cell_info(aquifer_index,2)-1 + lateral.PARENT.STATVAR_PRIVATE.overlap_info{i,1}(:,1);
                
                for n=1:size(flux_between_aquifers,1)
                    
                    flux_in(range_own_aquifer(n,1),1) = flux_in(range_own_aquifer(n,1),1) + double(flux_between_aquifers(n,1) > 0) .* flux_between_aquifers(n,1);
                    flux_out(range_own_aquifer(n,1),1) = flux_out(range_own_aquifer(n,1),1) - double(flux_between_aquifers(n,1) < 0) .* flux_between_aquifers(n,1);
                    
                    T_water = double(flux_between_aquifers(n,1) < 0) .* lateral.PARENT.STATVAR.T_water(range_own_aquifer(n,1),1);  %no advection of heat from resvervoir
                    flux_energy_in(range_own_aquifer(n,1),1) = flux_energy_in(range_own_aquifer(n,1),1) + double(flux_between_aquifers(n,1) > 0) .* flux_between_aquifers(n,1) .* ...
                        (double(T_water>0) .* lateral.PARENT.CONST.c_w + double(T_water<0) .* lateral.PARENT.CONST.c_i) .* T_water;
                    flux_energy_out(range_own_aquifer(n,1),1) = flux_energy_out(range_own_aquifer(n,1),1) - double(flux_between_aquifers(n,1) < 0) .* flux_between_aquifers(n,1) .* ...
                        (double(T_water>0) .* lateral.PARENT.CONST.c_w + double(T_water<0) .* lateral.PARENT.CONST.c_i) .* T_water;
                end

            end
            
            lateral.STATVAR.flux_in = flux_in;
            lateral.STATVAR.flux_out = flux_out;
            lateral.STATVAR.flux_energy_in = flux_energy_in;
            lateral.STATVAR.flux_energy_out = flux_energy_out;
            
            
            lateral.PARENT.STATVAR.water_flux = (flux_in - flux_out).* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
            lateral.PARENT.STATVAR.water_flux_energy = (flux_energy_in - flux_energy_out) .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
            

            
            
%             for i=1:size(enumeration_pressure_unknown,1)
%                 K_prime_times_head(enumeration_pressure_unknown(i,1),:) = K_prime(enumeration_pressure_unknown(i,1),:) .* unknown_heads(i,1);
%             end
%             
%             %flux_matrix = K_prime.*(head-head');
%             flux_matrix = K_prime_times_head - K_prime_times_head';
%             
%             lateral.STATVAR.flux_matrix = flux_matrix;
%             
%             total_fluxes = sum(flux_matrix,1)';
%             
%             fluxes = {};
%             offset = 0;
%             for i = 1:lateral.PARENT.PARA.num_realizations
%                 fluxes{i,1} = total_fluxes(offset+1:offset+number_of_aquifers_per_worker(i,1),1);
%                 offset = offset + number_of_aquifers_per_worker(i,1);
%             end
%             
%             lateral.STATVAR.fluxes = fluxes;
%             lateral.STATVAR.number_of_aquifers_per_worker = number_of_aquifers_per_worker;

            
            labBarrier;
        end

        
        
        function lateral = push(lateral, tile)
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
                    CURRENT = lateral3D_push_water_general_aquifer(CURRENT, lateral);
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
        
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE_LATERAL(ground);
            
            ground.PARA.class_category = 'LATERAL_IA';
            
            ground.PARA.options = [];
            ground.PARA.STATVAR = [];
            
            ground.PARA.default_value.hardBottom_cutoff = {0.03};
            ground.PARA.comment.hardBottom_cutoff = {'hard bottom  = no water flow if saturated and water content below [vol. water content, -]'};
            
            ground.PARA.default_value.ia_time_increment = {0.25};
            ground.PARA.comment.ia_time_increment ={'time step [days], must be multiple of of LATERAL class timestep'};
        end
        
    end
    
end


