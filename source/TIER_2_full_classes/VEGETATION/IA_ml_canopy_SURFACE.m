classdef IA_ml_canopy_SURFACE < IA_BASE

    
    properties
        SURFACE %points to the SURFACE class (GROUND or SNOW)
        VEGATION
    end
    
    methods

        function get_boundary_condition_u_vegetation_GROUND_surface(ia_vegetation_surface)
            ground = ia_vegetation_surface.SURFACE;
            vegetation = ia_vegetation_surface.VEGATION;
            forcing = vegetation.ForcingV;
            
            % fluxes from ground and snow overwritten by fluxes from vegetation module
            ground.STATVAR.Qh = vegetation.STATVAR.vegetation.mlcanopyinst.shsoi;
            ground.STATVAR.Qe = vegetation.STATVAR.vegetation.mlcanopyinst.lhsoi;
            ground.TEMP.F_ub = vegetation.STATVAR.vegetation.mlcanopyinst.gsoi .* ground.STATVAR.area(1,1);
            ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) + ground.TEMP.F_ub; %BE CAREFUL, THIS MAY DOUBLE-COUNT THE GROUND HEAT FLUX!
            
            
            vegetation.STATVAR.vegetation.soilvar.t_top_surfacecell = ground.STATVAR.T(1) + 273.15;
            vegetation.STATVAR.vegetation.soilvar.dz_topsurfacecell = ground.STATVAR.layerThick(1);
            vegetation.STATVAR.vegetation.soilvar.thk_topsurfacecell = ground.STATVAR.thermCond(1);

            
            
        end
        
        function get_boundary_condition_u_vegetation_SNOW_surface(ia_vegetation_surface)
            ground = ia_vegetation_surface.SURFACE;
            vegetation = ia_vegetation_surface.VEGATION;
            forcing = vegetation.ForcingV;
            

            
            
        end
        
    end
end

