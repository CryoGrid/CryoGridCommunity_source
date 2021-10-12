classdef DA_particle_batch_smoother < matlab.mixin.Copyable

    
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
        end
        
        function da = provide_CONST(da)
            
        end
        
        function da = provide_STATVAR(da)
            
        end
        
        function da = finalize_init(da, tile)
            da.TILE = tile;
            da.TEMP.index_next_obs = [];
            da.TEMP.time_next_obs = [];
            for i=1:size(da.PARA.observation_files,1)
                temp=load([da.PARA.observation_paths{i,1} da.PARA.observation_files{i,1}], 'OBS');
                da.STATVAR.obs_time{i,1} = temp.OBS.time;
                da.STATVAR.observations{i,1} = temp.OBS.observations;
                da.STATVAR.obs_error_variance{i,1} = temp.OBS.obs_error_variance;
                da.STATVAR.modeled_obs{i,1} = da.STATVAR.observations{i,1}.*NaN;
                
                da.OBS_OP{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.observable_classes{i,1}){da.PARA.observable_classes_index(i,1)});     
                da.TEMP.index_next_obs = [da.TEMP.index_next_obs; 1];
                da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(1,1)];
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
            
            
        end
        
        function da = DA_step(da, tile)
            if tile.t>= da.DA_TIME
                %loop over all observation data sets
                for i=1:size(da.STATVAR.obs_time,1)
                    if tile.t>= da.TEMP.time_next_obs(i,1)
                        da.STATVAR.modeled_obs{i,1}(da.TEMP.index_next_obs(i,1),1) = observable_operator(da.OBS_OP{i,1}, tile);
                        if da.TEMP.index_next_obs(i,1) < size(da.STATVAR.observations{i,1}, 1) %end of observations reached
                            da.TEMP.index_next_obs(i,1) = da.TEMP.index_next_obs(i,1) + 1;
                            da.TEMP.time_next_obs(i,1) = da.STATVAR.obs_time{i,1}(da.TEMP.index_next_obs(i,1),1);
                        else
                            da.TEMP.time_next_obs(i,1) = tile.FORCING.PARA.end_time
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
                    modeled_obs = [modeled_obs; da.STATVAR.modeled_obs{i,1}];
                end
                data_package = pack(da, data_package, 'modeled_obs', modeled_obs);
                da.ENSEMBLE.modeled_obs = repmat(modeled_obs, 1, da.TILE.RUN_INFO.PARA.number_of_tiles) .* NaN;
                da.ENSEMBLE.modeled_obs(:, da.TILE.RUN_INFO.PARA.worker_number) = modeled_obs;
                
                ensemble_param = []; %gather ensemble parameters in one vector
                for i=1:size(da.TILE.ENSEMBLE, 1)
                    variables = fieldnames(da.TILE.ENSEMBLE{i,1}.STATVAR);
                    for j=1:size(variables,1)
                        ensemble_param = [ensemble_param; da.TILE.ENSEMBLE{i,1}.STATVAR.(variables{j,1})];
                    end
                end
                data_package = pack(da, data_package, 'ensemble_param', ensemble_param);
                da.ENSEMBLE.ensemble_param = repmat(ensemble_param, 1, da.TILE.RUN_INFO.PARA.number_of_tiles) .* NaN;
                da.ENSEMBLE.ensemble_param(:, da.TILE.RUN_INFO.PARA.worker_number) = ensemble_param;
                
                %send
                for i = 1:da.TILE.RUN_INFO.PARA.number_of_tiles
                    if i~=da.TILE.RUN_INFO.PARA.worker_number
                        labSend(data_package, i, 1);
                    end
                end
                
                for i = 1:da.TILE.RUN_INFO.PARA.number_of_tiles
                    if i~=da.TILE.RUN_INFO.PARA.worker_number
                        data_package_in = labReceive(i, 1);
                        if ~isempty(data_package_in)
                            da = unpack(da, data_package_in, i); %read received column vector and transform into ENSEMBLE
                        end
                    end
                end
                
                %actual DA here
                da = PBS(da);
                
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
            No=size(da.STATVAR.observations{1,1},1);
            if numel(da.STATVAR.obs_error_variance{1,1})==No
                if size(da.STATVAR.obs_error_variance{1,1},2)==No
                    Rinv=da.STATVAR.obs_error_variance{1,1}.^(-1);
                else
                    Rinv=(da.STATVAR.obs_error_variance{1,1}').^(-1);
                end
            elseif numel(da.STATVAR.obs_error_variance{1,1})==1
                Rinv=(1/da.STATVAR.obs_error_variance{1,1}).*ones(1,No);
            else
                error('Expected numel(R)=No or scalar R')
            end
            
            % Calculate the likelihood.
            Inn=repmat(da.STATVAR.observations{1,1},1,size(da.ENSEMBLE.modeled_obs,2))-da.ENSEMBLE.modeled_obs;   % Innovation.
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

    end
end

