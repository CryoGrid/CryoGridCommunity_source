classdef GEOTHERMAL_Davies2013 < matlab.mixin.Copyable

    properties
        PARA
        CONST
        STATVAR
        TEMP
    end

    
    methods
        
        function geothermal = provide_PARA(geothermal)
            geothermal.PARA.geothermal_file = [];
            geothermal.PARA.geothermal_path = [];
        end
        
        function geothermal = provide_STATVAR(geothermal)

        end
        
        function geothermal = provide_CONST(geothermal)
            
        end
        
        function geothermal = finalize_init(geothermal)
            
        end
        

        function heat_flux = get_geothermal_heatflux(geothermal, run_info)
            load([geothermal.PARA.geothermal_path 'heatflowDavies2013.mat']);
            heat_flux = interp2(lat_mat, lon_mat, heat_flow_mat, run_info.STATVAR.latitude, run_info.STATVAR.longitude)./1000;  %map in W/mK
        end
        
        
    end
end

