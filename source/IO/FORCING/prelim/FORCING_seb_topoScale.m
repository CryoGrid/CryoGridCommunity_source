classdef FORCING_seb_topoScale < FORCING_base & READ_FORCING_mat & PROCESS_FORCING_topoScale
    
    properties
        
    end
    
    methods
        function forcing = provide_PARA(forcing)
            % INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.
            
            forcing.PARA.filename = [];   %filename of Matlab file containing forcing data
            forcing.PARA.forcing_path = [];
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.all_rain_T = [];     % Temperature above which all precipitation is considered as rain
            forcing.PARA.all_snow_T = [];     % Temperature below which all precipitation is considered as snow
            forcing.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.airT_height = [];  % height above ground at which air temperature (and wind speed!) from the forcing data are applied.
        end
        
        
        function forcing = provide_CONST(forcing)
            forcing.CONST.Tmfw = [];
            forcing.CONST.sigma = [];
        end
        
        
        function forcing = provide_STATVAR(forcing)
            
        end
        
        
        function forcing = finalize_init(forcing, tile)
            %add something like "if the ending if filename is mat, use
            %read_mad era4d, if it is nc, use a single nc-file, if it is a
            %directory, assume that it is piecewise nc-files, which must be
            %read by Kris' scripts -> in this case, set a variable as
            %switchso that the forcing reads a new data chunk and reapplies TopoScale as soon as it
            %is over the end of the chunk
            
            forcing.TEMP.era = read_mat_ERA4D(forcing, [forcing.PARA.forcing_path forcing.PARA.filename]);
            
            forcing = process_topoScale(forcing, tile);
            forcing.TEMP.era = []; 
            
            forcing = split_precip_Tair(forcing); % distinguish snow-/rainfall

            forcing.DATA.rainfall =  forcing.DATA.rainfall.*forcing.PARA.rain_fraction;
            forcing.DATA.snowfall =  forcing.DATA.snowfall.*forcing.PARA.snow_fraction;
            
            forcing = check_and_correct(forcing); % Remove known errors
            forcing = set_start_and_end_time(forcing); % assign start/end time
            forcing = initialize_TEMP(forcing);
            
            
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            forcing = interpolate_forcing@FORCING_base(forcing, tile);
            
            forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
            forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
        end
        
    end
    
end