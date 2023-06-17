classdef PERTURB_forcing_simple < matlab.mixin.Copyable

    
    properties
        PARA
        CONST
        STATVAR
    end
    
    methods
        function perturb = provide_PARA(perturb)
            perturb.PARA.Tair_bias = [];
            perturb.PARA.Sin_rel_error = [];
            perturb.PARA.snow_fraction = [];
            perturb.PARA.rain_fraction = [];
            perturb.PARA.all_rain_T = [];
            perturb.PARA.all_snow_T = [];
        end
        
        function perturb = provide_CONST(perturb)
            perturb.CONST.sigma = [];
        end
        
        function perturb = provide_STATVAR(perturb)
            
        end
        
        function perturb = finalize_init(perturb, tile)

        end
        
        function forcing = perturb_forcing(perturb, forcing)
            sky_emissivity = forcing.TEMP.Lin ./ (forcing.TEMP.Tair+273.15).^4 ./ perturb.CONST.sigma;
            forcing.TEMP.Tair = forcing.TEMP.Tair + perturb.PARA.Tair_bias;
            forcing.TEMP.Lin = sky_emissivity .* perturb.CONST.sigma .* (forcing.TEMP.Tair+273.15).^4;
            forcing.TEMP.Sin = forcing.TEMP.Sin .* (1 + perturb.PARA.Sin_rel_error);
            
            total_precip = forcing.TEMP.rainfall + forcing.TEMP.snowfall;
            forcing.TEMP.rainfall = total_precip .* (double(forcing.TEMP.Tair >= perturb.PARA.all_rain_T) + ...
                double(forcing.TEMP.Tair > perturb.PARA.all_snow_T & forcing.TEMP.Tair < perturb.PARA.all_rain_T) .* ...
                (forcing.TEMP.Tair - perturb.PARA.all_snow_T) ./ max(1e-12, (perturb.PARA.all_rain_T - perturb.PARA.all_snow_T)));
            forcing.TEMP.snowfall = total_precip - forcing.TEMP.rainfall;
            forcing.TEMP.snowfall = forcing.TEMP.snowfall .* perturb.PARA.snow_fraction;
            forcing.TEMP.rainfall = forcing.TEMP.rainfall .* perturb.PARA.rain_fraction;

        end

    end
end

