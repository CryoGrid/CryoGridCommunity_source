% %========================================================================
% CryoGrid TIER1 library class, functions related to vegetation
% R. B. Zweigel, June 2021
%========================================================================

classdef VEGETATION < BASE
    
    methods
        
        function canopy = get_heat_capacity_canopy_leaves(canopy)
            % As in Bonan et al. (2018), Modeling canopy-induced turbulence in the Earth system
            % Assumes stems have same heat capacity as leaves
            L = canopy.STATVAR.LAI;
            S = canopy.STATVAR.SAI;
            SLA = canopy.PARA.SLA;
            f_c = canopy.PARA.f_carbon;
            f_w = canopy.PARA.f_water;
            c_w = canopy.CONST.c_w./canopy.CONST.rho_w;
            c_dry = c_w ./ 3; 
            
            Ma = 1/SLA; % leaf dry mass per unit area
            c_leaf_areal = c_dry.*Ma./f_c + c_w.*(f_w./(1-f_w)).*Ma./f_c;
            canopy.STATVAR.c_canopy = (L+S).*c_leaf_areal.*canopy.STATVAR.area(1); % [J/K]
        end
        
        function canopy = get_heat_capacity_canopy(canopy)
            % As get_get_heat_capacity_canopy(..), but with trunk heat
            % capacity from Swenson et al. (2018)
            kv = canopy.PARA.kv;
            D_bh = canopy.PARA.D_bh;
            h_tree = canopy.STATVAR.layerThick;
            N_tree = canopy.PARA.N_tree;
            rho_wood = canopy.PARA.rho_wood;
            L = canopy.STATVAR.LAI;
            SLA = canopy.PARA.SLA;
            f_c = canopy.PARA.f_carbon;
            f_w = canopy.PARA.f_water;
            c_w = canopy.CONST.c_w./canopy.CONST.rho_w; % [J/(kg*K)]
            c_dry = c_w ./ 3; % specific heat capacity of dry biomass
            
            Ma = 1/SLA; % leaf dry mass per unit area
            c_leaf_areal = c_dry.*Ma./f_c + c_w.*(f_w./(1-f_w)).*Ma./f_c;
            
            V_tree = kv.*pi.*(D_bh/2).^2.*h_tree; % trunk volume
            M_tree =  N_tree.*rho_wood.*V_tree; % tree dry mass per area 
            c_stem = ( c_dry + c_w.*(f_w./(1-f_w)) ).*M_tree;
            
            canopy.STATVAR.c_canopy = (L.*c_leaf_areal + c_stem).*canopy.STATVAR.area(1); % [J/K]
        end
        
        function canopy = get_T_simpleVegetatation(canopy)
            % Disregards water in canopy
            canopy.STATVAR.T = canopy.STATVAR.energy ./ canopy.STATVAR.c_canopy;
        end
        
        function canopy = get_E_simpleVegetation(canopy)
            canopy.STATVAR.energy = canopy.STATVAR.T .* canopy.STATVAR.c_canopy; % [J]
        end
        
        function canopy = get_E_water_vegetation(canopy)
            canopy.STATVAR.water = double(canopy.STATVAR.T>=0).*canopy.STATVAR.waterIce; % [m3]
            canopy.STATVAR.ice = double(canopy.STATVAR.T<0).*canopy.STATVAR.waterIce; % [m3]
            c_canopy = canopy.STATVAR.c_canopy; % [J/K]
            cp_waterIce = (canopy.STATVAR.water.*canopy.CONST.c_w + canopy.STATVAR.ice.*canopy.CONST.c_i); % [J/K]
            
            canopy.STATVAR.energy = canopy.STATVAR.T.*(c_canopy + cp_waterIce) - double(canopy.STATVAR.T<0).*canopy.STATVAR.waterIce.*canopy.CONST.L_f; % [J]
        end
        
        function canopy = get_T_water_vegetation(canopy)
            Lf = canopy.CONST.L_f; % [J/m3]
            c_i = canopy.CONST.c_i; % [J/m3/K]
            c_w = canopy.CONST.c_w; % [J/m3/K]
            c_canopy = canopy.STATVAR.c_canopy; % [J/K]
            L = canopy.STATVAR.LAI;
            S = canopy.STATVAR.SAI;
            W_max = canopy.PARA.Wmax; 
            
            e_frozen = -Lf.*canopy.STATVAR.waterIce;
            
            canopy.STATVAR.T = double(canopy.STATVAR.energy < e_frozen).*(canopy.STATVAR.energy - e_frozen)./(c_i.*canopy.STATVAR.waterIce + c_canopy) ...
                + double(canopy.STATVAR.energy > 0).* canopy.STATVAR.energy./(c_w.*canopy.STATVAR.waterIce + c_canopy);
            canopy.STATVAR.ice = double(canopy.STATVAR.energy <= e_frozen).*canopy.STATVAR.waterIce + double(canopy.STATVAR.energy > e_frozen & canopy.STATVAR.energy < 0) .* canopy.STATVAR.energy./(-Lf);
            canopy.STATVAR.water = double(canopy.STATVAR.energy >= 0).*canopy.STATVAR.waterIce + double(canopy.STATVAR.energy > e_frozen & canopy.STATVAR.energy < 0).*(canopy.STATVAR.energy - e_frozen)./Lf;
            
            canopy.STATVAR.f_wet = ( canopy.STATVAR.waterIce ./ (W_max.*canopy.STATVAR.area(1).*(L+S)) ).^(2/3);
            canopy.STATVAR.f_wet = min(1,canopy.STATVAR.f_wet);
            canopy.STATVAR.f_dry = (1-canopy.STATVAR.f_wet).*L./(L+S);
        end
        
        function timestep = get_timestep_canopy_T(canopy)
            d_energy = canopy.TEMP.d_energy;
            c_canopy =  canopy.STATVAR.c_canopy; % [J/m3]
            cp_waterIce = (canopy.STATVAR.water.*canopy.CONST.c_w + canopy.STATVAR.ice.*canopy.CONST.c_i); % [J/K]
            
            timestep = canopy.PARA.dT_max ./ ( abs(d_energy)./(c_canopy + cp_waterIce) );
        end
        
        function canopy = get_z0_d_vegetation(canopy) % roughness length of vegetated surface (CLM5)
            L = canopy.STATVAR.LAI; % Leaf area index
            S = canopy.PARA.SAI; % Stem area index
            z0g = get_z0_surface(canopy.NEXT); % Roughness lenght of ground/snow surface
            R_z0 = canopy.PARA.R_z0; % Ratio of momentum roughness length to canopy height
            R_d = canopy.PARA.R_d; % Ratio of displacement height to canopy height
            z_top = sum(canopy.STATVAR.layerThick); % canopy height
            
            V = ( 1-exp(-1.*min(L+S,2))) ./ (1-exp(-2)); % Eq. 5.127
            z0 = exp( V.*log(z_top.*R_z0) + (1-V).*log(z0g) ); % Eq. 5.125
            d = z_top.*R_d.*V; % Eq. 5.126
            canopy.STATVAR.z0 = z0;
            canopy.STATVAR.d = d;
        end
        
        function canopy = add_canopy(canopy)
            canopy.STATVAR.LAI = canopy.PARA.LAI;
            canopy.STATVAR.emissivity = 1 - exp(-canopy.STATVAR.LAI-canopy.STATVAR.SAI);
            canopy = get_heat_capacity_canopy(canopy);
            canopy = get_E_water_vegetation(canopy); % derive energy from temperature
            canopy = get_z0_d_vegetation(canopy);
        end
        
        function canopy = remove_canopy(canopy)
            canopy.STATVAR.LAI = 0;
            canopy.STATVAR.emissivity = 1 - exp(-canopy.STATVAR.LAI-canopy.STATVAR.SAI);
            canopy = get_heat_capacity_canopy(canopy);
            canopy = get_E_water_vegetation(canopy); % derive energy from temperature            
            canopy = get_z0_d_vegetation(canopy);
        end
        
%         function canopy = distribute_roots(canopy)
%             beta = canopy.PARA.beta_root;
%             dz = canopy.GROUND.STATVAR.layerThick;
%             z = cumsum(dz);
%            
%             % Root fraction per soil layer
%             f_root = beta.^([0; z(1:end-1)].*100) - beta.^(z*100); % Eq. 11.1
%             
%             canopy.GROUND.STATVAR.f_root = f_root;
%         end
%         
%         function stresses = get_soil_moisture_stress(canopy)
%             psi = canopy.GROUND.STATVAR.waterPotential;
%             f_root = canopy.GROUND.STATVAR.f_root;
% %             layerThick = ia_seb_water.NEXT.STATVAR.layerThick;
%             psi_wilt = canopy.PARA.psi_wilt;
%             
%             stresses = sum(f_root.*min(1,(psi./psi_wilt)));
%         end
    end
end

