%========================================================================
% CryoGrid TIER1 INTERACTION (IA) class for functions related to surface
% energy balance extending through the top class (e.g. canopy/vegetation)
% R. B. Zweigel, August 2021
%========================================================================

classdef IA_SEB < IA_BASE
    
    methods
        
        function ia_seb = get_boundary_condition_Qh_m(ia_seb, tile)
            forcing = tile.FORCING;
            stratigraphy1 = ia_seb.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb.NEXT; % ground
            
            cp = stratigraphy2.CONST.cp;
            kappa = stratigraphy2.CONST.kappa;
            
            uz = forcing.TEMP.wind;
            z =  stratigraphy2.PARA.airT_height;
            z0 = stratigraphy2.PARA.z0;
            Tz = forcing.TEMP.Tair;
            TForcing = stratigraphy2.STATVAR.T(1);
            Lstar = stratigraphy2.STATVAR.Lstar;
            p = forcing.TEMP.p;
            
            Tz=Tz+forcing.CONST.Tmfw;
            TForcing=TForcing+forcing.CONST.Tmfw;
            rho = rho_air(stratigraphy2,p, Tz);
            
            Q_h  = -rho.*cp.*kappa.* uz.*kappa./(log(z./z0)- psi_M(stratigraphy2, z./Lstar, z0./Lstar)) .* (Tz-TForcing)./(log(z./z0)- psi_H(stratigraphy2, z./Lstar, z0./Lstar));
            stratigraphy2.STATVAR.Qh = Q_h;
        end
        
        function ia_seb = get_boundary_condition_Qe_m(ia_seb, tile)
            forcing = tile.FORCING;
            stratigraphy1 = ia_seb.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb.NEXT; % ground
            
            stratigraphy2 = Q_evap_CLM4_5(stratigraphy2, forcing);
            
            stratigraphy2.TEMP.d_water_ET(1,1) = stratigraphy2.TEMP.d_water_ET(1,1) -  stratigraphy2.STATVAR.evap.* stratigraphy2.STATVAR.area(1,1); %in m3 water per sec, put everything in uppermost grid cell
            stratigraphy2.TEMP.d_water_ET_energy(1,1) = stratigraphy2.TEMP.d_water_ET_energy(1,1) -  stratigraphy2.STATVAR.evap_energy.* stratigraphy2.STATVAR.area(1,1);
        end
        
    end
end
