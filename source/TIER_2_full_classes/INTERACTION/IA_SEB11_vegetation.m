%========================================================================
% CryoGrid INTERACTION (IA) class for surface energy balance of a GROUND
% class below a shading VEGETATION class
% R. B. Zwegiel, August 2021
%========================================================================

classdef IA_HEAT01_SEB11_vegetation < IA_SEB % IA_WATER & IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat_seb, tile)
            % Get sensible- & latent fluxes as if lower class was top_class
            ia_heat_seb = get_boundary_condition_Qh_m(ia_heat_seb, tile);
            ia_heat_seb = get_boundary_condition_Qe_m(ia_heat_seb, tile);
            
            % add fluxes to uppermost cell (ratiative fluxes are added by penetration)
            ia_heat_seb.NEXT.TEMP.d_energy(1) = ia_heat_seb.NEXT.TEMP.d_energy(1) + (-ia_heat_seb.NEXT.STATVAR.Qh - ia_heat_seb.NEXT.STATVAR.Qe) .* ia_heat_seb.NEXT.STATVAR.area(1);
        end
        

        
    end
end