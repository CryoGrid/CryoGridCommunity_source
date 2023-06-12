classdef GEOTHERMAL_Davies2013 < matlab.mixin.Copyable

    properties
        PARENT
        CONST
        STATVAR
        PARA
    end

    
    methods
        
        function heat_flux = provide_PARA(heat_flux)
            heat_flux.PARA.geothermal_file = [];
            heat_flux.PARA.geothermal_path = [];
        end
        
        function heat_flux = provide_STATVAR(heat_flux)

        end
        
        function heat_flux = provide_CONST(heat_flux)
            
        end
        
        function heat_flux = finalize_init(heat_flux)
            
        end
        
        function heat_flux = load_data(heat_flux)

          %  heat_flux.PARENT.STATVAR.geothermal = heat_flux.PARENT.STATVAR.latitude .* 0+0.05;
            
            load([heat_flux.PARA.geothermal_path heat_flux.PARA.geothermal_file]);
            heat_flux.PARENT.STATVAR.geothermal = interp2(lat_mat, lon_mat, heat_flow_mat, heat_flux.PARENT.STATVAR.latitude, heat_flux.PARENT.STATVAR.longitude)./1000;  %map in W/mK
        end
        
        
    end
end

