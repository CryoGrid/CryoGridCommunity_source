classdef DA_particle_filter < matlab.mixin.Copyable

    
    properties
        TILE
        OBS_OP
        PARA
        CONST
        STATVAR
        TEMP
        DA_TIME        
    end
    
    methods
        function da = provide_PARA(da)
            da.PARA.observation_files = [];
            da.PARA.observation_paths = [];
            da.PARA.observable_classes = [];
            da.PARA.observable_classes_index = []; %must all have the same length, i.e. each observational data set requires one observable class 
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
                da.STATVAR.modeled_obs{i,1} = da.STATVAR.observations{i,1}.*NaN;
                
                da.OBS_OP{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.observable_classes{i,1}){da.PARA.observable_classes_index(i,1)});     
                da.TEMP.index_next_obs = [da.TEMP.index_next_obs; 1];
                da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(1,1)];
            end

            da.DA_TIME = min(da.TEMP.time_next_obs);
        end
        
        function da = DA_step(da, tile)
            if tile.t>= da.DA_TIME
                %loop over all observation data sets
                for i=1:size(da.STATVAR.obs_time,1)
                    if tile.t>= da.TEMP.time_next_obs(i,1)
                        da.STATVAR.modeled_obs{i,1}(da.TEMP.index_next_obs(i,1),1) = observable_operator(da.OBS_OP{i,1}, tile);
                        if da.TEMP.index_next_obs(i,1) < size(da.STATVAR.observations{i,1}, 1)
                            da.TEMP.index_next_obs(i,1) = da.TEMP.index_next_obs(i,1) + 1;
                            da.TEMP.time_next_obs(i,1) = da.STATVAR.obs_time{i,1}(da.TEMP.index_next_obs(i,1),1);
                        else
                            da.TEMP.time_next_obs(i,1) = tile.FORCING.PARA.end_time
                        end
                    end
                end
                
                da.DA_TIME = min(da.TEMP.time_next_obs);
                
                %add sync and DA step here
                %...
                
            end
        end

    end
end

