%========================================================================
% CryoGrid INTERACTION (IA) class for interaction between a vegetaion class
% and an belowlying soil class. Only for subsurface - canopy interactions,
% such as transpiration, no surface energy/water balance.
% R. B. Zwegiel, February 2022
%========================================================================

classdef IA_VEGETATION_CLM5 < IA_WATER
    
    properties
    end
    
    methods
        
        function ia_soil = distribute_roots_CLM5(ia_soil, tile)
            % Vertical root distribution as in 2.11.1.1 in CLM5 documentation
            beta = ia_soil.PREVIOUS.PARA.beta_root; % Root distribution parameter
            dz = ia_soil.NEXT.STATVAR.layerThick;
            z = cumsum(dz);
           
            % Root fraction per soil layer
            f_root = beta.^([0; z(1:end-1)].*100) - beta.^(z*100); % Eq. 11.1
            
            ia_soil.NEXT.STATVAR.f_root = f_root;
        end
        
        function beta = get_soil_moisture_stress(ia_soil, tile)
            % as described by Sinclair (2005), and Verhoef & Egea (2014)
            psi = ia_soil.NEXT.STATVAR.waterPotential;
            f_root = ia_soil.NEXT.STATVAR.f_root;
            psi_wilt = ia_soil.PREVIOUS.PARA.psi_wilt;
            
            beta = sum(f_root.*max(0,1-psi./psi_wilt));
        end
       
        function ia_soil = get_water_transpiration(ia_soil)
            ia_soil = get_water_transpiration@IA_WATER(ia_soil);
        end
        
    end
end