classdef IA_ml_canopy_GROUND < IA_BASE

    
    properties
        GROUND
        VEGETATION
    end
    
    methods
        function get_transpiration(ia_vegetation_ground)
            
            ground = ia_vegetation_surface.SURFACE;
            vegetation = ia_vegetation_surface.VEGATION;
            
            vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
            
            
            ground.TEMP.d_water(1,1) = ground.TEMP.d_water(1,1) - vegetation.STATVAR.vegetation.soilvar.transp_per_layer(1,1) .* ground.STATVAR.area(1,1);
            ground.TEMP.d_water(2:3,1) = ground.TEMP.d_water(2:3,1) - vegetation.STATVAR.vegetation.soilvar.transp_per_layer(2,1)/2 .* ground.STATVAR.area(2:3,1);
            ground.TEMP.d_water(4:9,1) = ground.TEMP.d_water(4:9,1) - vegetation.STATVAR.vegetation.soilvar.transp_per_layer(3,1)/6 .* ground.STATVAR.area(4:9,1);

        end
        
        function get_evaporation(ia_vegetation_ground)
            
            ground = ia_vegetation_surface.SURFACE;
            vegetation = ia_vegetation_surface.VEGATION;

            vegetation.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* vegetation.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
            
            ground.TEMP.d_water(1,1) = ground.TEMP.d_water(1,1) - vegetation.STATVAR.vegetation.mlcanopyinst.etsoi .* ground.STATVAR.area(1,1);

        end
    end
end

