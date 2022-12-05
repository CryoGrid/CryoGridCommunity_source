% base class build a model tile

classdef RUN_3D_PARALLEL < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        TILE
        %pprovider
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)
            
            run_info.PARA.number_of_cores = [];
            run_info.PARA.number_of_tiles = []; %3;
            run_info.PARA.param_file_number = []; %[1;2;3];
            run_info.PARA.run_mode = 'parallel';   % or 'sequential'

            run_info.PARA.connected = [];
            run_info.PARA.contact_length = [];
            run_info.PARA.distance = [];

            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];

        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        
%         function run_info = initialize_excel(run_info)
%             
%         end
        
        function run_info = finalize_init(run_info)
            if isempty(run_info.PARA.run_mode) || isnan(run_info.PARA.run_mode)
                run_info.PARA.run_mode = 'parallel';
            end
        end
        
        
        function [run_info, tile] = run_model(run_info, run_flag)
            %this could first open spmd and assign run_number depending on
            %worker, then do another round of pprovider
            %it could also do a loop over different tile representing
            %different sections of the run, e.g. initial inial init, spin-up, actual run 
            %
            % The run_info and tile instances returned are not the ones
            % actually used in the runs, it is the template instances
            % that were passed to the function.
            %
            % run_flag is a flag to indicate whether the model should be
            % run or only initialized. Setting it to false is not very
            % meaningful in this class, since initialization will be lost
            % due to the parallelization. But in sequential mode, the last
            % tile wil be returned in initialized state.
            
            if ~exist('run_flag', 'var')
                % if run_flag is not passed, default to true
                run_flag = true;
            end

            err_out = cell(run_info.PARA.number_of_tiles);
            
            tStart = datetime(datestr(now));

            if strcmpi(run_info.PARA.run_mode, 'parallel')
                this_pool = gcp('nocreate'); 
                if isempty(this_pool)
                    this_pool = parpool([1 run_info.PARA.number_of_cores]);
                    disp(['Using new parpool with ' num2str(this_pool.NumWorkers) ' workers.'])
                else
                    disp(['Using existing parpool with ' num2str(this_pool.NumWorkers) ' workers.'])
                end
                    
                
                parfor tile_id = 1:run_info.PARA.number_of_tiles
                    % make copy of template run_info, to modify in this process
                    this_run_info = copy(run_info);
                        this_run_info.PARA.worker_number = tile_id; % assign id
                    
                    % prepare error collection in case of exceptions
                    err_out{tile_id}.tile_id = tile_id;
                    err_out{tile_id}.OK = true;
    
                    try
                        % initialize and run tile instance
                        [out_run_info, out_tile] = kernel_run_model(this_run_info, run_flag);
                    catch ME
                        % catch and store error for saving after pool completes
                        error_timestamp = now;
    
                        err_out{tile_id}.tile_id = tile_id;
                        err_out{tile_id}.OK = false;
                        err_out{tile_id}.run_info = this_run_info;
                        err_out{tile_id}.timestamp = error_timestamp;
                        err_out{tile_id}.MException = ME;
    
                        % we cannot use save inside parfor, so we have to store
                        % the information and save it after.
                    end
                end
                tile = run_info.TILE;

            elseif strcmpi(run_info.PARA.run_mode, 'sequential')
                for tile_id = 1:run_info.PARA.number_of_tiles
                    this_run_info = copy(run_info);
                    this_run_info.PARA.worker_number = tile_id; % assign id
                    
                    % prepare error collection in case of exceptions
                    err_out{tile_id}.tile_id = tile_id;
                    err_out{tile_id}.OK = true;
                    
                    try
                        % initialize and run tile instance
                        [out_run_info, out_tile] = kernel_run_model(this_run_info, run_flag);
                    catch ME
                        % catch and store error for saving after pool completes
                        error_timestamp = now;
    
                        err_out{tile_id}.tile_id = tile_id;
                        err_out{tile_id}.OK = false;
                        err_out{tile_id}.run_info = this_run_info;
                        err_out{tile_id}.timestamp = error_timestamp;
                        err_out{tile_id}.MException = ME;
                    end
                end
                
                run_info = out_run_info;
                tile = out_tile;

            end

            deltatime = datetime(datestr(now))-tStart;
            fprintf('Elapsed time: ');
            disp(deltatime);

            % now save any error logs
            for tid = 1:length(run_info.PARA.number_of_tiles)
                if ~err_out{tid}.OK
                    error_log_file = [err_out{tid}.run_info.PPROVIDER.PARA.result_path, ...
                                      'errorlog_', ...
                                      datestr(err_out{tid}.timestamp, 'yyyymmdd_HHMMSS'), '__', ...
                                      err_out{tid}.run_info.PPROVIDER.PARA.run_name, ...
                                      '__tile_', num2str(err_out{tid}.tile_id), ...
                                      '.mat'];
                    info_out = err_out{tid};
                    save(error_log_file, '-struct', 'info_out');
                end
            end                

        end
 
        
        function [run_info, tile] = setup_run(run_info)
            %update the worker-specific name of the parameter file and the run
            if ~isfield(run_info.PARA, 'worker_number')
                run_info.PARA.worker_number = 1;
            end
            
            run_info.PPROVIDER = update_parameter_file(run_info.PPROVIDER, run_info.PARA.param_file_number(run_info.PARA.worker_number,1));
            run_info.PPROVIDER = update_run_name(run_info.PPROVIDER, run_info.PARA.worker_number);
                            
            %read worker-specific parameter file
            run_info.PPROVIDER = read_parameters(run_info.PPROVIDER);
                            
            run_info = customize(run_info);
            
            %tile = run_info.PPROVIDER.FUNCTIONAL_CLASSES.TILE{run_info.PARA.tile_number,1};
            tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});

            tile.RUN_INFO = run_info;
            run_info.TILE = tile;

            tile = finalize_init(tile);
        end


        function [run_info, tile] = kernel_run_model(run_info, run_flag)
            % This is the actual normal run_model method. It is extracted
            % in separate method to more easily enclose its execution in
            % a try-catch block in the new run_model method.

            [run_info, tile] = setup_run(run_info);

            if run_flag
                tile = run_model(tile);  %time integration
            end
        end
    

        function run_info = customize(run_info)
            %FUNCTION TO BE EDITED BY USER - when inheriting from this class
            %here customizations can be done by directly writing pprovider
            %one can for example change parameters in the different
            %subsurface classes
            %run_info.PPROVIDER.FUNCTIONAL_CLASSES.XX = YY;
            %
            %use run_info.PARA.worker_number for customization

        end
        
    end
end



