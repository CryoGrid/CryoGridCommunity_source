classdef FORCING_slope_seb_surfaceLevel < FORCING_base & READ_FORCING_mat & PROCESS_FORCING_topoScale
    
    properties
        
    end
    
    methods
        function forcing = provide_PARA(forcing)
            
            forcing.PARA.filename = [];   %filename of Matlab file containing forcing data
            forcing.PARA.forcing_path = [];
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.all_rain_T = [];     % Temperature above which all precipitation is considered as rain
            forcing.PARA.all_snow_T = [];     % Temperature below which all precipitation is considered as snow
            forcing.PARA.albedo_surrounding_terrain = [];
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
            
            forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));
            
            forcing.TEMP.era = read_mat_ERA4D(forcing, [forcing.PARA.forcing_path forcing.PARA.filename]);
            
            forcing = interpolate_sl(forcing, tile);
            forcing.TEMP.era = []; 
            
            forcing = split_precip_Tair(forcing); % distinguish snow-/rainfall

            forcing.DATA.rainfall =  forcing.DATA.rainfall.*forcing.PARA.rain_fraction;
            forcing.DATA.snowfall =  forcing.DATA.snowfall.*forcing.PARA.snow_fraction;
            
            forcing = check_and_correct(forcing); % Remove known errors
%             forcing = set_start_and_end_time(forcing); % assign start/end time
            forcing = initialize_TEMP(forcing);
            forcing = initialize_TEMP_slope(forcing);

            forcing = reduce_precip_slope(forcing, tile);
            
            forcing = SolarAzEl(forcing, tile);
            %make own function?
            if ~isfield(forcing.DATA, 'S_TOA')
                mu0=max(sind(forcing.DATA.sunElevation),0); % Trunacte negative values.
                sunset=mu0<cosd(89);%(mu0==0); % Sunset switch.
                forcing.DATA.S_TOA = 1370.*mu0;
            end
            
            forcing = split_Sin(forcing);
            forcing = terrain_corr_Sin_dif(forcing, tile);
            forcing = reproject_Sin_dir(forcing, tile);
            forcing = terrain_shade(forcing, tile);
            forcing = terrain_corr_Lin(forcing, tile);
            forcing.DATA.Sin = forcing.DATA.Sin_dir + forcing.DATA.Sin_dif;
            
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            forcing = interpolate_forcing@FORCING_base(forcing, tile);
            
            forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
            forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
              
        end
        
    end
    
end