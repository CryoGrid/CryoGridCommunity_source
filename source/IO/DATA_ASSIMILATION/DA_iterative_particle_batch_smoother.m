classdef DA_iterative_particle_batch_smoother < matlab.mixin.Copyable

    
    properties
        TILE
        OBS_OP
        PARA
        CONST
        STATVAR
        TEMP
        DA_TIME     
        DA_STEP_TIME
        ENSEMBLE
    end
    
    methods
        function da = provide_PARA(da)
            da.PARA.observation_files = [];
            da.PARA.observation_paths = [];
            da.PARA.observable_classes = [];
            da.PARA.observable_classes_index = []; %must all have the same length, i.e. each observational data set requires one observable class 
            da.PARA.assimilation_frequency = []; %year, month or day
            da.PARA.assimilation_interval = []; %number of years, months or days
            da.PARA.assimilation_date = []; %specific date in case of years or months
            da.PARA.start_assimilation_period = []; %Hlist, date from when the assimilation is started, i.e. the initial state
            da.PARA.ensemble_variables = [];
            da.PARA.learning_coefficient_iteration = [];
            da.PARA.learning_coefficient_forward = [];
            da.PARA.max_iterations = [];
        end
        
        function da = provide_CONST(da)
            
        end
        
        function da = provide_STATVAR(da)
            
        end
        
        function da = finalize_init(da, tile)
            da.TILE = tile;
%             offset = length(num2str(da.TILE.PARA.worker_number))+1;
            da.TEMP.run_name = da.TILE.PARA.run_name;%(1, 1:length(da.TILE.PARA.run_name)-offset); %original results folder w.o worker number
            da.TEMP.first_obs_index =[];
            da.TEMP.index_next_obs = [];
            da.TEMP.time_next_obs = [];
            for i=1:size(da.PARA.observation_files,1)
                temp=load([da.PARA.observation_paths{i,1} da.PARA.observation_files{i,1}], 'OBS');
                da.STATVAR.obs_time{i,1} = temp.OBS.time;
                da.STATVAR.observations{i,1} = temp.OBS.observations;
                da.STATVAR.obs_variance{i,1} = temp.OBS.obs_variance;
                da.STATVAR.modeled_obs{i,1} = da.STATVAR.observations{i,1}.*NaN;
                da.ENSEMBLE.weights = repmat(1./da.TILE.PARA.ensemble_size, 1, da.TILE.PARA.ensemble_size); %set equal weights before 1st DA step
                
                da.OBS_OP{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.observable_classes{i,1}){da.PARA.observable_classes_index(i,1)});     
                da.TEMP.first_obs_index = [da.TEMP.first_obs_index; find(da.STATVAR.obs_time{i,1} > tile.t, 1)];
                da.TEMP.index_next_obs = [da.TEMP.index_next_obs; da.TEMP.first_obs_index(i,1)]; %start with 
                da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(da.TEMP.first_obs_index(i,1),1)];
            end

            da.DA_TIME = min(da.TEMP.time_next_obs);
            
            if strcmp(da.PARA.assimilation_frequency, 'year')
                current_year = str2num(datestr(tile.t, 'yyyy'));
                da.DA_STEP_TIME = datenum([da.PARA.assimilation_date num2str(current_year + da.PARA.assimilation_interval)], 'dd.mm.yyyy');
            elseif strcmp(da.PARA.assimilation_frequency, 'month')
                current_year = str2num(datestr(tile.t, 'yyyy'));
                current_month = str2num(datestr(tile.t, 'mm'));
                da.DA_STEP_TIME = datenum(current_year, current_month + da.PARA.assimilation_interval, da.PARA.assimilation_date);
            elseif strcmp(da.PARA.assimilation_frequency, 'day')
                da.DA_STEP_TIME = tile.t + da.PARA.assimilation_interval;
            end
          
            if isempty(da.PARA.start_assimilation_period) || sum(isnan(da.PARA.start_assimilation_period))>0
                da = save_state(da, tile);
                da.TEMP.last_assimilation_date = tile.t+1; %start very close to initial state
                da.TEMP.assimilation_started = 0;
            else
                da.TEMP.assimilation_started = 0;
                da.TEMP.last_assimilation_date = datenum(da.PARA.start_assimilation_period(1,1), da.PARA.start_assimilation_period(2,1), da.PARA.start_assimilation_period(3,1));
            end
            
            da.TEMP.num_iterations = 1;
        end
        
        function da = DA_step(da, tile)
            %save the state and the ensemble variables when the DA begins,
            %this step is redone every time the DA is happy and moves on in time
            if ~da.TEMP.assimilation_started && tile.t>=da.TEMP.last_assimilation_date
                da.TEMP.last_assimilation_date = tile.t;
                da.TEMP.assimilation_started = 1;
                da.ENSEMBLE.weights_old = da.ENSEMBLE.weights;
                da = save_state(da, tile); %save states at the start of the assimilation period 
                for i=1:size(da.PARA.ensemble_variables,1)
                    da.TEMP.old_ensemble.(da.PARA.ensemble_variables{i,1}) = get_variable_info(tile.ENSEMBLE, da.PARA.ensemble_variables{i,1});
                end
            end
            
            if da.TEMP.assimilation_started && tile.t>= da.DA_TIME
                %loop over all observation data sets
                for i=1:size(da.STATVAR.obs_time,1)
                    if tile.t>= da.TEMP.time_next_obs(i,1)
                        da.STATVAR.modeled_obs{i,1}(da.TEMP.index_next_obs(i,1),1) = observable_operator(da.OBS_OP{i,1}, tile);
                        disp(da.STATVAR.modeled_obs{i,1}(da.TEMP.index_next_obs(i,1),1) )
                        if da.TEMP.index_next_obs(i,1) < size(da.STATVAR.observations{i,1}, 1) 
                            da.TEMP.index_next_obs(i,1) = da.TEMP.index_next_obs(i,1) + 1;
                            da.TEMP.time_next_obs(i,1) = da.STATVAR.obs_time{i,1}(da.TEMP.index_next_obs(i,1),1);
                        else %end of observations reached
                            da.TEMP.time_next_obs(i,1) = tile.FORCING.PARA.end_time +1;
                        end
                    end
                    
                end
                
                da.DA_TIME = min(da.TEMP.time_next_obs);
                
            end
            
            if tile.t>=da.DA_STEP_TIME
                labBarrier;
                %synchronize
                data_package = [];
                
                modeled_obs = [];%gather modeled observations in one vector
                for i=1:size(da.STATVAR.modeled_obs,1) 
                    modeled_obs = [modeled_obs; da.STATVAR.modeled_obs{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,1)];  %ONLY USE THE PART IN THE OBS INTERVAL, ALSO MAKE A SIMILAR VECTOR FOR OBSERVATIONS AND VARIANCES
                end
                data_package = pack(da, data_package, 'modeled_obs', modeled_obs);
                da.ENSEMBLE.modeled_obs = repmat(modeled_obs, 1, da.TILE.PARA.ensemble_size) .* NaN;
                da.ENSEMBLE.modeled_obs(:, da.TILE.PARA.worker_number) = modeled_obs;
                
                ensemble_param = []; %gather ensemble parameters in one vector
                variables = fieldnames(da.TILE.ENSEMBLE.STATVAR);
                for j=1:size(variables,1)
                    ensemble_param = [ensemble_param; da.TILE.ENSEMBLE.STATVAR.(variables{j,1})];
                end
                data_package = pack(da, data_package, 'ensemble_param', ensemble_param);
                da.ENSEMBLE.ensemble_param = repmat(ensemble_param, 1, da.TILE.PARA.ensemble_size) .* NaN;
                da.ENSEMBLE.ensemble_param(:, da.TILE.PARA.worker_number) = ensemble_param;
                
                %send
                for i = 1:da.TILE.PARA.ensemble_size
                    if i~=da.TILE.PARA.worker_number
                        labSend(data_package, i, 1);
                    end
                end
                
                for i = 1:da.TILE.PARA.ensemble_size
                    if i~=da.TILE.PARA.worker_number
                        data_package_in = labReceive(i, 1);
                        if ~isempty(data_package_in)
                            da = unpack(da, data_package_in, i); %read received column vector and transform into ENSEMBLE
                        end
                    end
                end
                
                %actual DA 
                da = PBS(da);
                
                %MAKE THIS PART OF the OUT class
                if da.TILE.PARA.worker_number==1
                    test= copy(da);
                    test.TILE = [];
                    save('test.mat', 'test')
                end
                
                %at this point the weights are known 
                %1. determine effective sample size
                %2. if effective sample size > threshold, save the
                %surviving ensemble members, read them into the ones w.o.
                %weight and then recompute ensemble parameters according to
                %defined protocol (learning or non-learning)
                %3. if effective ensemble size < threshold, do not save
                %new, but reset timestamps, again read the surviving
                %ensemble members from the last round, and recompute the
                %ensemble parameters based on the weights of the surviving ensemble members 
                
                effective_ensemble_size = 1./sum(da.ENSEMBLE.weights.^2);
                
                if effective_ensemble_size >= da.TILE.PARA.ensemble_size/5 || da.TEMP.num_iterations>=da.PARA.max_iterations
                    
                    da.TEMP.num_iterations = 1;
                    
                    %resampling from weights
                    rng(tile.t+25) %use current time as seed for randum number generator, i.e. same sequence of random numbers will be generated by each worker
                    [weights_sorted, posi] = sort(da.ENSEMBLE.weights);
                    weights_sorted_cum = cumsum(weights_sorted);
                    weights_sorted_cum = [0 weights_sorted_cum(1, 1:end-1)];
                    resample_posi=sum(double(repmat(rand(da.TILE.PARA.ensemble_size,1),1,da.TILE.PARA.ensemble_size) > repmat(weights_sorted_cum, da.TILE.PARA.ensemble_size,1)),2);
                    get_from_new_worker = [];
                    for i=1:da.TILE.PARA.ensemble_size
                        get_from_new_worker=[get_from_new_worker; posi(resample_posi(i))];
                    end
                    
                    %save the stratigraphy vector and all the other TILE info
                    %in file
                    if sum(get_from_new_worker==da.TILE.PARA.worker_number)>0
                        da = save_state(da, tile);
                    end
                    
                    labBarrier;
                    
                    %read the new stratigraphy and info from file
                    temp=load([tile.PARA.result_path  da.TEMP.run_name '/tile_' num2str(get_from_new_worker(da.TILE.PARA.worker_number,1)) '.mat']);
                    variables = fieldnames(temp.state);
                    for i=1:size(variables,1)
                        if ~isempty(temp.state.(variables{i,1}))
                            tile.(variables{i,1}) = temp.state.(variables{i,1});
                        end
                    end
                    
                    labBarrier;
                    if sum(get_from_new_worker==da.TILE.PARA.worker_number)>0
                        delete([tile.PARA.result_path da.TEMP.run_name '/tile_' num2str(da.TILE.PARA.worker_number) '.mat']);
                    end
                    
                    
                    %recalculate ensemble parameters
                    ENSEMBLE_variables = fieldnames(tile.ENSEMBLE.STATVAR);
                    for i=1:size(ENSEMBLE_variables)
                        if any(strcmp(ENSEMBLE_variables{i,1}, da.PARA.ensemble_variables)) 
                            %find out if that variable is to be optimized by this DA class and check learning coefficient if this ensemble member is supposed to learn
                            if any(strcmp(tile.ENSEMBLE.PARA.gaussian_variables_name, ENSEMBLE_variables{i,1})) %gaussian variable
                                pos = find(strcmp(tile.ENSEMBLE.PARA.gaussian_variables_name, ENSEMBLE_variables{i,1}));
                                if da.PARA.learning_coefficient_forward >= rand()
                                    tile.ENSEMBLE.PARA.gaussian_variables_center(pos,1) = sum(da.ENSEMBLE.ensemble_param(i,:) .* da.ENSEMBLE.weights);
                                    if effective_ensemble_size >= 3
                                        a = round(da.ENSEMBLE.weights .* da.TILE.PARA.ensemble_size.*100);
                                        b=[];
                                        for k=1:length(a)
                                            b=[b repmat(da.ENSEMBLE.ensemble_param(i,k), 1, a(1,k))];
                                        end
                                        tile.ENSEMBLE.PARA.gaussian_variables_width(pos,1) = std(b);
                                    end
                                else
                                    tile.ENSEMBLE.PARA.gaussian_variables_center(pos,1) = da.TEMP.old_ensemble.(ENSEMBLE_variables{i,1}).gaussian_variables_center;
                                    tile.ENSEMBLE.PARA.gaussian_variables_width(pos,1) = da.TEMP.old_ensemble.(ENSEMBLE_variables{i,1}).gaussian_variables_width;
                                end
                            end
                            if any(strcmp(tile.ENSEMBLE.PARA.boxcar_variables_name, ENSEMBLE_variables{i,1})) %boxcar variable
                                pos = find(strcmp(tile.ENSEMBLE.PARA.boxcar_variables_name, ENSEMBLE_variables{i,1}));
                                if da.PARA.learning_coefficient_forward >= rand()
                                    if effective_ensemble_size >= 3
                                        mean_value = sum(da.ENSEMBLE.ensemble_param(i,:) .* da.ENSEMBLE.weights);
                                        range = tile.ENSEMBLE.PARA.boxcar_variables_upper_bound(pos,1) - tile.ENSEMBLE.PARA.boxcar_variables_lower_bound(pos,1);
                                        tile.ENSEMBLE.PARA.boxcar_variables_lower_bound(pos,1) = mean_value - range/2;
                                        tile.ENSEMBLE.PARA.boxcar_variables_upper_bound(pos,1) = mean_value + range/2;
                                    else
                                        valid_pos = find(da.ENSEMBLE.weights(1,:) > 1./(100.*da.TILE.PARA.ensemble_size));
                                        mini = min(da.ENSEMBLE.ensemble_param(i,valid_pos));
                                        maxi = max(da.ENSEMBLE.ensemble_param(i,valid_pos));
                                        range = maxi-mini;
                                        tile.ENSEMBLE.PARA.boxcar_variables_lower_bound(pos,1) = mini - range/2;
                                        tile.ENSEMBLE.PARA.boxcar_variables_upper_bound(pos,1) = maxi + range/2;
                                    end
                                else
                                    tile.ENSEMBLE.PARA.boxcar_variables_lower_bound(pos,1) = da.TEMP.old_ensemble.(ENSEMBLE_variables{i,1}).boxcar_variables_lower_bound;
                                    tile.ENSEMBLE.PARA.boxcar_variables_upper_bound(pos,1) = da.TEMP.old_ensemble.(ENSEMBLE_variables{i,1}).boxcar_variables_upper_bound;
                                end
                            end
                        end
                    end
                    da.TILE.ENSEMBLE = recalculate_ensemble_parameters_after_DA(da.TILE.ENSEMBLE, tile, da.PARA.ensemble_variables);
                    
                    %assign next DA_STEP_TIME
                    if strcmp(da.PARA.assimilation_frequency, 'year')
                        current_year = str2num(datestr(da.DA_STEP_TIME, 'yyyy'));
                        da.DA_STEP_TIME = datenum([da.PARA.assimilation_date num2str(current_year + da.PARA.assimilation_interval)], 'dd.mm.yyyy');
                    elseif strcmp(da.PARA.assimilation_frequency, 'month')
                        current_year = str2num(datestr(da.DA_STEP_TIME, 'yyyy'));
                        current_month = str2num(datestr(da.DA_STEP_TIME, 'mm'));
                        da.DA_STEP_TIME = datenum(current_year, current_month + da.PARA.assimilation_interval, da.PARA.assimilation_date);
                    elseif strcmp(da.PARA.assimilation_frequency, 'day')
                        da.DA_STEP_TIME = da.DA_STEP_TIME + da.PARA.assimilation_interval;
                    end
                    
                    da.TEMP.first_obs_index =[];
                    da.TEMP.index_next_obs = [];
                    da.TEMP.time_next_obs = [];
                    for i=1:size(da.PARA.observation_files,1)
                        if ~isempty(find(da.STATVAR.obs_time{i,1} > tile.t, 1))
                            da.TEMP.first_obs_index = [da.TEMP.first_obs_index; find(da.STATVAR.obs_time{i,1} > tile.t, 1)];
                            da.TEMP.index_next_obs = [da.TEMP.index_next_obs; da.TEMP.first_obs_index(i,1)]; %start with
                            da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(da.TEMP.first_obs_index(i,1),1)];
                        else
                            da.TEMP.first_obs_index = [da.TEMP.first_obs_index; size(da.STATVAR.observations{i,1}, 1)];
                            da.TEMP.index_next_obs = [da.TEMP.index_next_obs; size(da.STATVAR.observations{i,1}, 1)];
                            da.TEMP.time_next_obs = [da.TEMP.time_next_obs; tile.FORCING.PARA.end_time + 1];
                        end
                    end
                    
                    da.DA_TIME = min(da.TEMP.time_next_obs);
                    da.TEMP.last_assimilation_date = tile.t;
                    da.ENSEMBLE.weights_old = da.ENSEMBLE.weights;
                    da = save_state(da, tile); %save new states at the start of the new assimilation period, so that it can be read again if ensemble is degenerate
                    for i=1:size(da.PARA.ensemble_variables,1)
                        da.TEMP.old_ensemble.(da.PARA.ensemble_variables{i,1}) = get_variable_info(tile.ENSEMBLE, da.PARA.ensemble_variables{i,1});
                    end
                    %assimilation successful, ensemble is not degenerate,
                    %move on in time
                else
                    %ensemble is degenerate, start over with old states 
                    disp('ensemble degenerate, one more time')
                    
                    da.TEMP.num_iterations = da.TEMP.num_iterations + 1;
                    
                    temp=load([tile.PARA.result_path da.TEMP.run_name '/tile_' num2str(da.TILE.PARA.worker_number) '.mat']);
                    variables = fieldnames(temp.state);
                    for i=1:size(variables,1)
                        if ~isempty(temp.state.(variables{i,1}))
                            tile.(variables{i,1}) = temp.state.(variables{i,1});
                        end
                    end
                    
                    %recalculate ensemble parameters, same as before, but
                    %using learning_coefficient_iteration
                    ENSEMBLE_variables = fieldnames(tile.ENSEMBLE.STATVAR);
                    for i=1:size(ENSEMBLE_variables)
                        if any(strcmp(ENSEMBLE_variables{i,1}, da.PARA.ensemble_variables)) 
                            %find out if that variable is to be optimized by this DA class and check learning coefficient if this ensemble member is supposed to learn
                            if any(strcmp(tile.ENSEMBLE.PARA.gaussian_variables_name, ENSEMBLE_variables{i,1})) %gaussian variable
                                pos = find(strcmp(tile.ENSEMBLE.PARA.gaussian_variables_name, ENSEMBLE_variables{i,1}));
                                if da.PARA.learning_coefficient_iteration >= rand() %carry the information from the DA step into the future 
                                    tile.ENSEMBLE.PARA.gaussian_variables_center(pos,1) = sum(da.ENSEMBLE.ensemble_param(i,:) .* da.ENSEMBLE.weights);
                                    if effective_ensemble_size >= 3
                                        a = round(da.ENSEMBLE.weights .* da.TILE.PARA.ensemble_size.*100);
                                        b=[];
                                        for k=1:length(a)
                                            b=[b repmat(da.ENSEMBLE.ensemble_param(i,k), 1, a(1,k))];
                                        end
                                        tile.ENSEMBLE.PARA.gaussian_variables_width(pos,1) = std(b);
                                    else
                                        tile.ENSEMBLE.PARA.gaussian_variables_width(pos,1) = tile.ENSEMBLE.PARA.gaussian_variables_width(pos,1)/2; %make reduction factor a PARA 
                                    end
                                else %go back to the previous values 
                                    tile.ENSEMBLE.PARA.gaussian_variables_center(pos,1) = da.TEMP.old_ensemble.(ENSEMBLE_variables{i,1}).gaussian_variables_center;
                                    tile.ENSEMBLE.PARA.gaussian_variables_width(pos,1) = da.TEMP.old_ensemble.(ENSEMBLE_variables{i,1}).gaussian_variables_width;
                                end
                            end
                            if any(strcmp(tile.ENSEMBLE.PARA.boxcar_variables_name, ENSEMBLE_variables{i,1})) %boxcar variable
                                pos = find(strcmp(tile.ENSEMBLE.PARA.boxcar_variables_name, ENSEMBLE_variables{i,1}));
                                if da.PARA.learning_coefficient_iteration >= rand() %carry the information from the DA step into the future
                                    if effective_ensemble_size <= 3
                                        mean_value = sum(da.ENSEMBLE.ensemble_param(i,:) .* da.ENSEMBLE.weights);
                                        range = tile.ENSEMBLE.PARA.boxcar_variables_upper_bound(pos,1) - tile.ENSEMBLE.PARA.boxcar_variables_lower_bound(pos,1);
                                        tile.ENSEMBLE.PARA.boxcar_variables_lower_bound(pos,1) = mean_value - range/4; %make reduction factor a PARA 
                                        tile.ENSEMBLE.PARA.boxcar_variables_upper_bound(pos,1) = mean_value + range/4;
                                    else
                                        valid_pos = find(da.ENSEMBLE.weights(1,:) > 1./(100.*da.TILE.PARA.ensemble_size));
                                        mini = min(da.ENSEMBLE.ensemble_param(i,valid_pos));
                                        maxi = max(da.ENSEMBLE.ensemble_param(i,valid_pos));
                                        range = maxi-mini;
                                        tile.ENSEMBLE.PARA.boxcar_variables_lower_bound(pos,1) = mini - range/4;
                                        tile.ENSEMBLE.PARA.boxcar_variables_upper_bound(pos,1) = maxi + range/4;  %make reduction factor a PARA 
                                    end
                                else
                                    tile.ENSEMBLE.PARA.boxcar_variables_lower_bound(pos,1) = da.TEMP.old_ensemble.(ENSEMBLE_variables{i,1}).boxcar_variables_lower_bound;
                                    tile.ENSEMBLE.PARA.boxcar_variables_upper_bound(pos,1) = da.TEMP.old_ensemble.(ENSEMBLE_variables{i,1}).boxcar_variables_upper_bound;
                                end
                            end
                        end
                    end
                    da.TILE.ENSEMBLE = recalculate_ensemble_parameters_after_DA(da.TILE.ENSEMBLE, tile, da.PARA.ensemble_variables);
                    
                    %reset timestamps, no need to reset timestamps in the
                    %CG stratigraphy since the old states are read in
                    tile.t = da.TEMP.last_assimilation_date;

                    da.TEMP.first_obs_index =[];
                    da.TEMP.index_next_obs = [];
                    da.TEMP.time_next_obs = [];
                    for i=1:size(da.PARA.observation_files,1)
                        if ~isempty(find(da.STATVAR.obs_time{i,1} > tile.t, 1))
                            da.TEMP.first_obs_index = [da.TEMP.first_obs_index; find(da.STATVAR.obs_time{i,1} > tile.t, 1)];
                            da.TEMP.index_next_obs = [da.TEMP.index_next_obs; da.TEMP.first_obs_index(i,1)]; %start with
                            da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(da.TEMP.first_obs_index(i,1),1)];
                        else
                            da.TEMP.first_obs_index = [da.TEMP.first_obs_index; size(da.STATVAR.observations{i,1}, 1)];
                            da.TEMP.index_next_obs = [da.TEMP.index_next_obs; size(da.STATVAR.observations{i,1}, 1)];
                            da.TEMP.time_next_obs = [da.TEMP.time_next_obs; tile.FORCING.PARA.end_time + 1];
                        end
                    end
                    
                    da.DA_TIME = min(da.TEMP.time_next_obs);
                    tile.OUT = reset_timestamp_out(tile.OUT,tile);
                    
                end
            end
            
        end
        
        function da = PBS(da)
            % Efficient implementation of the Particle Batch Smoother
            % presented in Margulis et al. (2015; JHM).
            % N.B. The observation errors are assumed to be uncorrelated (diagonal R)
            % and Gaussian.
            %
            % Dimensions: No = Number of observations in the batch to assimilate.
            %             Np = Number of parameters to update.
            %             Ne = Number of ensemble members (particles).
            %
            % -----------------------------------------------------------------------
            % Inputs:
            %
            %
            % HX   => No x Ne matrix containing an ensemble of Ne predicted
            %         observation column vectors each with No entries.
            %
            % Y     => No x 1 vector containing the batch of (unperturbed) observations.
            %
            % R     => No x No observation error variance matrix; this may also be
            %         specified as a scalar corresponding to the constant variance of
            %         all the observations in the case that these are all from the same
            %         instrument.
            %
            % -----------------------------------------------------------------------
            % Outputs:
            %
            % w     => 1 x Ne vector containing the ensemble of posterior weights,
            %         the prior weights are implicitly 1/N_e.
            %
            % -----------------------------------------------------------------------
            % See e.g. https://jblevins.org/log/log-sum-exp for underflow issue.
            %
            % Code by Kristoffer Aalstad (Feb. 2019)
            
            
            
            % Calculate the diagonal of the inverse obs. error covariance.
            observations = [];
            obs_variance = [];
            for i=1:size(da.STATVAR.observations,1)
                observations = [observations; da.STATVAR.observations{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,1)];
                obs_variance = [obs_variance; da.STATVAR.obs_variance{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,1)];
            end
            da.ENSEMBLE.observations = observations;
            da.ENSEMBLE.obs_variance = obs_variance;
            
            No=size(observations,1);
            Rinv=(obs_variance').^(-1);
%             if numel(da.STATVAR.obs_error_variance{1,1})==No
%                 if size(da.STATVAR.obs_error_variance{1,1},2)==No
%                     Rinv=da.STATVAR.obs_error_variance{1,1}.^(-1);
%                 else
%                     Rinv=(da.STATVAR.obs_error_variance{1,1}').^(-1);
%                 end
%             elseif numel(da.STATVAR.obs_error_variance{1,1})==1
%                 Rinv=(1/da.STATVAR.obs_error_variance{1,1}).*ones(1,No);
%             else
%                 error('Expected numel(R)=No or scalar R')
%             end
            
            % Calculate the likelihood.
            Inn=repmat(observations,1,size(da.ENSEMBLE.modeled_obs,2))-da.ENSEMBLE.modeled_obs;   % Innovation.
            EObj=Rinv*(Inn.^2);                     % [1 x Ne] ensemble objective function.
            LLH=-0.5.*EObj; % log-likelihoods.
            normc=logsumexp(da, LLH,2);
            
            
            % NB! The likelihood coefficient (1/sqrt(2*pi...)) is
            % omitted because it drops out in the normalization
            % of the likelihood. Including it (very small term) would lead
            % to problems with FP division.
            
            
            % Calculate the posterior weights as the normalized likelihood.
            logw = LLH-normc;
            da.ENSEMBLE.weights = exp(logw); % Posterior weights.
            
            % Need "log-sum-exp" trick to overcome numerical issues for small R/large
            % number of obs.
            
        end
            
            
            
            
            
        function s = logsumexp(da, a, dim)
            % Returns log(sum(exp(a),dim)) while avoiding numerical underflow.
            % Default is dim = 1 (columns).
            % logsumexp(a, 2) will sum across rows instead of columns.
            % Unlike matlab's "sum", it will not switch the summing direction
            % if you provide a row vector.
            
            % Written by Tom Minka
            % (c) Microsoft Corporation. All rights reserved.
            
            if nargin < 2
                dim = 1;
            end
            
            % subtract the largest in each column
            y = max(a,[],dim);
            dims = ones(1,ndims(a));
            dims(dim) = size(a,dim);
            a = a - repmat(y, dims);
            s = y + log(sum(exp(a),dim));
            i = find(~isfinite(y));
            if ~isempty(i)
                s(i) = y(i);
            end
        end
        
        
        function data_package = pack(da, data_package, var_name, var) %transform into column vector ready to send
                %variables{i,1}
                data_package=[data_package; size(var_name,2); double(var_name)']; % # of characters followed by characters as doubles
                data_package=[data_package; size(var,1); var]; % # of entries followed by values
        end
        
        function da = unpack(da, data_package, received_from_worker) %read received column vector and transform into STATVAR
            i=1;
            while i<=size(data_package,1)
               variable_name = char(data_package(i+1:i+data_package(i,1),1)');
               i = i + data_package(i,1) + 1;
               da.ENSEMBLE.(variable_name)(:,received_from_worker) = data_package(i+1:i+data_package(i,1),1);
               i = i + data_package(i,1) + 1;
            end
        end
        
        function da = save_state(da, tile)
            state = copy(tile);
            variables = fieldnames(state);
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'LATERAL') && ~strcmp(variables{i,1}, 'TOP') && ~strcmp(variables{i,1}, 'BOTTOM') && ~strcmp(variables{i,1}, 'TOP_CLASS') && ~strcmp(variables{i,1}, 'BOTTOM_CLASS')
                    state.(variables{i,1}) = [];
                end
            end
            
            save([tile.PARA.result_path da.TEMP.run_name '/tile_' num2str(da.TILE.PARA.worker_number) '.mat'], 'state');

        end

    end
end

