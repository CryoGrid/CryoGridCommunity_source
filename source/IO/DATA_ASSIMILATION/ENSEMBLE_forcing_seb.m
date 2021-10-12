classdef ENSEMBLE_forcing_seb < matlab.mixin.Copyable

    
    properties
        TILE
        PARA
        CONST
        STATVAR
    end
    
    methods
        function ensemble = provide_PARA(ensemble)
            ensemble.PARA.absolute_error_Tair = [];
            ensemble.PARA.min_snow_fraction = [];
            ensemble.PARA.max_snow_fraction = [];
            ensemble.PARA.relative_error_Sin = [];     
            ensemble.PARA.all_snow_T = [];
            ensemble.PARA.all_rain_T = [];
        end
        
        function ensemble = provide_CONST(ensemble)
            ensemble.CONST.sigma = [];
        end
        
        function ensemble = provide_STATVAR(ensemble)
            
        end
        
        function ensemble = finalize_init(ensemble, tile)
            ensemble.TILE = tile;
            
            %assign constant values throughout entire run
            ensemble.STATVAR.Tair_bias = randn(tile.RUN_INFO.PARA.number_of_tiles, 1) ;
            ensemble.STATVAR.Tair_bias = ensemble.STATVAR.Tair_bias(tile.RUN_INFO.PARA.worker_number, 1) .* ensemble.PARA.absolute_error_Tair;
            ensemble.STATVAR.Sin_rel_error = randn(tile.RUN_INFO.PARA.number_of_tiles, 1);
            ensemble.STATVAR.Sin_rel_error = ensemble.STATVAR.Sin_rel_error(tile.RUN_INFO.PARA.worker_number, 1) .* ensemble.PARA.relative_error_Sin;
            ensemble.STATVAR.snow_fraction = rand(tile.RUN_INFO.PARA.number_of_tiles, 1);
            ensemble.STATVAR.snow_fraction = ensemble.PARA.min_snow_fraction + ensemble.STATVAR.snow_fraction(tile.RUN_INFO.PARA.worker_number, 1) .* (ensemble.PARA.max_snow_fraction - ensemble.PARA.min_snow_fraction);
            
        end
        
        function ensemble = ensemble_step(ensemble, tile)
            sky_emissivity = tile.FORCING.TEMP.Lin ./ (tile.FORCING.TEMP.Tair+273.15).^4 ./ ensemble.CONST.sigma;
            tile.FORCING.TEMP.Tair = tile.FORCING.TEMP.Tair + ensemble.STATVAR.Tair_bias;
            tile.FORCING.TEMP.Lin = sky_emissivity .* ensemble.CONST.sigma .* (tile.FORCING.TEMP.Tair+273.15).^4;
            tile.FORCING.TEMP.Sin = tile.FORCING.TEMP.Sin .* ensemble.STATVAR.Sin_rel_error;
            
            total_precip = tile.FORCING.TEMP.rainfall + tile.FORCING.TEMP.snowfall;
            tile.FORCING.TEMP.rainfall = total_precip .* (double(tile.FORCING.TEMP.Tair >= ensemble.PARA.all_rain_T)  + ...
                double(tile.FORCING.TEMP.Tair > ensemble.PARA.all_snow_T & tile.FORCING.TEMP.Tair < ensemble.PARA.all_rain_T) .* ...
                (1 - (tile.FORCING.TEMP.Tair - ensemble.PARA.all_snow_T) ./ max(1e-12, (ensemble.PARA.all_rain_T - ensemble.PARA.all_snow_T))));
            tile.FORCING.TEMP.snowfall = total_precip .* (double(tile.FORCING.TEMP.Tair <= ensemble.PARA.all_snow_T)  + ...
                double(tile.FORCING.TEMP.Tair > ensemble.PARA.all_snow_T & tile.FORCING.TEMP.Tair <= ensemble.PARA.all_rain_T) .* ...
                (tile.FORCING.TEMP.Tair - ensemble.PARA.all_snow_T) ./ max(1e-12, (ensemble.PARA.all_rain_T - ensemble.PARA.all_snow_T)));
            tile.FORCING.TEMP.snowfall = tile.FORCING.TEMP.snowfall .* ensemble.STATVAR.snow_fraction;

        end
        
        
        
    end
end

