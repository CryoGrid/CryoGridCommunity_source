
classdef TILE_ESA_CCI_preproc < matlab.mixin.Copyable
    
    properties
        PARA
        RUN_INFO
        BUILDER
        FORCING
        CONST
        OUT        
        PREPROC_CLASS
        
        
        t        
        timestep
        
    end
    
    
    methods
        
       

        function tile = provide_PARA(tile)
                        
            tile.PARA.preproc_class =[];
            tile.PARA.preproc_class_index =[];
            tile.PARA.forcing_class = [];
            tile.PARA.forcing_class_index = [];

            tile.PARA.out_class = [];
            tile.PARA.out_class_index = [];
            
            tile.PARA.time_interval = []; %in days
            
        end
        
        function tile = provide_CONST(tile)

        end
        
        function tile = provide_STATVAR(tile)

        end

        
        %assemble the stratigraphy
        function tile = finalize_init(tile)
            tile.RUN_INFO.TILE = tile;
            
            %2. forcing -> special forcing class required
            tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
            tile.FORCING = finalize_init(tile.FORCING, tile);
            
            tile.PREPROC_CLASS = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.preproc_class){tile.PARA.preproc_class_index,1});
            tile.PREPROC_CLASS = finalize_init(tile.PREPROC_CLASS, tile); %must contain function that attached MODIS LST class from RUN_INFO
            
            %10. assign time, etc.
            tile.t = tile.FORCING.PARA.start_time;
            %set fixed timestep!
            tile.timestep = tile.PREPROC_CLASS.PARA.timestep;

 
            %12. assign OUT classes
            tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
            tile.OUT = finalize_init(tile.OUT, tile);
            


        end
        

        function tile = interpolate_forcing_tile(tile)
             tile.FORCING = interpolate_forcing(tile.FORCING, tile);
        end
        

        
        function tile = store_OUT_tile(tile)
            tile.OUT = store_OUT(tile.OUT, tile);
        end        
        
        
        
        function tile = run_model(tile)
            
            %=========================================================================
            %TIME INTEGRATION
            %=========================================================================
            while tile.t < tile.FORCING.PARA.end_time
                
                
                
                CURRENT = tile.PREPROC_CLASS;
                
                disp(datestr(tile.t))
                
                %interpolate focing data to time t
                tile = interpolate_forcing_tile(tile);  %get_ERA_downscaled
                
                %upper boundar condition (uppermost class only)
                CURRENT = get_boundary_condition_u(CURRENT, tile); %get_MODIS
                
%                 %lower boundary condition (lowermost class)
%                 CURRENT = get_boundary_condition_l(CURRENT,  tile);
                
                %calculate spatial derivatives
                CURRENT = get_derivatives_prognostic(CURRENT, tile);  %calculate snowfall and melt

%                 %prognostic step - integrate prognostic variables in time
%                 CURRENT = advance_prognostic(CURRENT, tile);  
                
                %diagnostic step - compute diagnostic variables
                CURRENT = compute_diagnostic(CURRENT, tile);  %merge MODIS and ERA

%                 %triggers
%                 CURRENT = check_trigger(CURRENT, tile);
                
                %update timestep
                old_year = str2num(datestr(tile.t, 'yyyy'));
                %update time variable t
                tile.t = tile.t + tile.timestep;
                new_year = str2num(datestr(tile.t, 'yyyy'));
                
                if new_year>old_year %reset to Jan 1 to start over
                    tile.t = datenum(new_year, 1, 1);
                end
                
                
                %model output
                tile = store_OUT_tile(tile);
            end
            
        end
        

    end
end



