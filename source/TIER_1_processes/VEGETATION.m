%========================================================================
% CryoGrid TIER1 library class, functions related to vegetation
% R. B. Zweigel, June 2021
%========================================================================

classdef VEGETATION < BASE
    
    methods
               
        function canopy = get_T_simpleVegetatation(canopy)
            % Supersimplified, assumes heat capacity is ONLY dependent on
            % LAI, disregards phase change etc.
            lai = canopy.PARA.LAI; % [m2/m2]
            leaf_cp_areal = canopy.PARA.leaf_cp_areal; % [J/m2/K]
            % cp = 1999; % value from Blanken et al. (1997) (as in Bonan et. (2018)) for aspen leaves, sla = 111 g/m2 and waterfraction of 80%
            
            cp = lai .* leaf_cp_areal; % [J/K]
            canopy.STATVAR.T = canopy.STATVAR.energy ./ cp;
        end
        
        function canopy = get_E_simpleVegetation(canopy)
            lai = canopy.PARA.LAI; % [m2/m2]
            leaf_cp_areal = canopy.PARA.leaf_cp_areal; % [J/m2/K]
            
            cp = lai .* leaf_cp_areal;
            canopy.STATVAR.energy = canopy.STATVAR.T .* cp;
        end
        
        function timestep = get_timestep_canopy_T(canopy)
            d_energy = canopy.TEMP.d_energy;
            lai = canopy.PARA.LAI; % [m2/m2]
            leaf_cp_areal = canopy.PARA.leaf_cp_areal; % [J/m2/K]
            
            cp = lai .* leaf_cp_areal;
            
            timestep = canopy.PARA.dT_max ./ (abs(d_energy)./cp);
        end
    end
end

