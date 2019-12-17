classdef IA_HEAT_VEGETATION < matlab.mixin.Copyable
    
    %> Simone: Interaction class between Vegetation and Ground
    %> Soil heat flux (W/m^2) up and down
    
    properties
        PREVIOUS
        NEXT
    end
    
    methods
        
        %         Heat:
        %           Vegetation -> Ground
        %               Ground heat flux
        %               Throughfall rain and snow (maybe serperate)
        %           Ground -> Vegetation
        %               T0 an Vegetation zusammen mit conductivity and capacity
        
        %         Moisture:
        %           Vegetation -> Ground
        %               Layer bis max Wurzeltiefe, für diese Layer müssen die Variablen aus dem Boden an die Vegetation gegeben werden.
        %               Wasser kommt nur aus Boden nicht von Vegetation in Boden zurück
        %               Die Variablen in Vegetation: watsat, wksat, bsw and smp_l (mm) sollten über IA class vom Boden an die Vegetation zurückgegeben werden.
        
        
        function  get_boundary_condition_m(ia_heat_vegetation)
            stratigraphy1 = ia_heat_vegetation.PREVIOUS;
            stratigraphy2 = ia_heat_vegetation.NEXT;
            
            % Flux aus canopy in Boden
            stratigraphy1.TEMP.F_lb = -stratigraphy1.STATVAR.vegetation.mlcanopyinst.gsoi;
            stratigraphy2.TEMP.F_ub = -stratigraphy1.TEMP.F_lb;
            
            % Flux aus Boden in Canopy (wahrscheinlich nicht nötig)
            stratigraphy1.STATVAR.vegetation.soilvar.thk = stratigraphy2.STATVAR.thermCond(1); % thermal conductivity
            stratigraphy1.STATVAR.vegetation.soilvar.cv = stratigraphy2.STATVAR.mineral(1) .* stratigraphy2.CONST.c_m + stratigraphy2.STATVAR.organic(1) .* stratigraphy2.CONST.c_o + stratigraphy2.STATVAR.water(1) .* stratigraphy2.CONST.c_w + stratigraphy2.STATVAR.ice(1) .* stratigraphy2.CONST.c_i; % heat capacity  wasser / Luft / ice
            stratigraphy1.STATVAR.vegetation.soilvar.t_soisno = stratigraphy2.STATVAR.T(1)+273.15;
            
            % Influx of water into the ground
            stratigraphy2.STATVAR.water(1) = stratigraphy1.STATVAR.vegetation.mlcanopyinst.qflx_prec_grnd_rain + stratigraphy1.STATVAR.vegetation.mlcanopyinst.qflx_prec_grnd_snow;
            
            % Ground water content / saturation / saturation vapor pressure --> Available water for plants
            % At the moment the root fraction per layer is 0.2 -> 5 layers
            
%             stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_vol = stratigraphy2.STATVAR.water(1);
%             stratigraphy1.STATVAR.vegetation.soilvar.watsat = stratigraphy2.STATVAR.water(1);
%             stratigraphy1.STATVAR.vegetation.soilvar.bsw = stratigraphy2.STATVAR.water(1);
%             stratigraphy1.STATVAR.vegetation.soilvar.smp_l = stratigraphy2.STATVAR.water(1); % [mm] % Soil layer matric potential (mm)
            
            ground.STATVAR.vegetation.soilvar.hksat        = 0.176; %3; %0.00001; %0.01; %3.;   %soilstate_inst.hksat_col      ;  % Soil layer hydraulic conductivity at saturation (mm H2O/s)
            ground.STATVAR.vegetation.soilvar.bsw          = 4.05;   %soilstate_inst.bsw_col        ;  % Soil layer Clapp and Hornberger "b" parameter
            ground.STATVAR.vegetation.soilvar.rootfr       = 1;  %soilstate_inst.rootfr_patch   ;  % Fraction of roots in each soil layer
            ground.STATVAR.vegetation.soilvar.h2osoi_vol   = stratigraphy2.STATVAR.water(1);  %waterstate_inst.h2osoi_vol_col;  % Soil layer volumetric water content (m3/m3)
            ground.STATVAR.vegetation.soilvar.h2osoi_ice   = stratigraphy2.STATVAR.ice(1);   %waterstate_inst.h2osoi_ice_col;  % Soil layer ice lens (kg/m2)
            
            ground.STATVAR.vegetation.soilvar.z = 0.05;     % Soil layer depth (m)
            ground.STATVAR.vegetation.soilvar.zi = 0.05;     % Soil layer depth at layer interface (m)
            
            ground.STATVAR.vegetation.soilvar.nsoi = 1;
            
            % Soil layer thickness (m)
            ground.STATVAR.vegetation.soilvar.dz = 0.05;
            
            % Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface)
            ground.STATVAR.vegetation.soilvar.z_plus_onehalf(1) = -ground.STATVAR.vegetation.soilvar.dz(1);
            for i = 2:ground.STATVAR.vegetation.soilvar.nsoi
                ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i) = ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i-1) - ground.STATVAR.vegetation.soilvar.dz(i);
            end
            
            % Soil depth (m) at center of layer i (negative distance from surface)
            
            ground.STATVAR.vegetation.soilvar.z(1) = 0.5 * ground.STATVAR.vegetation.soilvar.z_plus_onehalf(1);
            for i = 2:ground.STATVAR.vegetation.soilvar.nsoi
                ground.STATVAR.vegetation.soilvar.z(i) = 0.5 * (ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i-1) + ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i));
            end
            %
            % %Thickness between between z(i) and z(i+1)
            
            for i = 1:ground.STATVAR.vegetation.soilvar.nsoi-1
                ground.STATVAR.vegetation.soilvar.dz_plus_onehalf(i) = ground.STATVAR.vegetation.soilvar.z(i) - ground.STATVAR.vegetation.soilvar.z(i+1);
            end
            ground.STATVAR.vegetation.soilvar.dz_plus_onehalf(ground.STATVAR.vegetation.soilvar.nsoi) = 0.5 * ground.STATVAR.vegetation.soilvar.dz(ground.STATVAR.vegetation.soilvar.nsoi);
            
            % --- Initial conditions
            
            % Initial soil temperature (K) and unfrozen and frozen water (kg H2O/m2)
            
            for i = 1:ground.STATVAR.vegetation.soilvar.nsoi
                
                % Temperature
                tmean = 15.0; %15.0;             % Mean daily air temperature (C)
                ground.STATVAR.vegetation.soilvar.tsoi(i) = tmean + 273.15;
                
                % Soil water at saturation (kg H2O/m2)
                
                % Watsat and soil_texture taken from sp_05_01.m
                % Volumetric soil water content at saturation (porosity)
                % (Clapp and Hornberger. 1978. Water Resources Research 14:601-604)
                
                ground.STATVAR.vegetation.soilvar.watsat = 0.3; %, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482];
                
                h2osoi_sat = ground.STATVAR.vegetation.soilvar.watsat(i) * 1000 * ground.STATVAR.vegetation.soilvar.dz(i);
                 
                % ground.STATVAR.vegetation.soilvar.soil_water_matric_potential = -(1/alpha)*(((0w - 0wr) / (0ws - 0wr))^(n /(1-n))-1)^(1/n);
                
                %soil parameters sandy clay loam (p.120, bonan book)
                ground.STATVAR.vegetation.soilvar.alpha = 5.9; %cm -> m  0.059;
                ground.STATVAR.vegetation.soilvar.n = 1.48;
                ground.STATVAR.vegetation.soilvar.Owr(i) = ground.STATVAR.vegetation.soilvar.watsat(i);
                ground.STATVAR.vegetation.soilvar.Ow(i) = ground.STATVAR.vegetation.soilvar.h2osoi_vol(i);
                ground.STATVAR.vegetation.soilvar.Ows = 0;
                
                ground.STATVAR.vegetation.soilvar.soil_water_matric_potential(i) = -(1/ground.STATVAR.vegetation.soilvar.alpha)*(((ground.STATVAR.vegetation.soilvar.Ow - ground.STATVAR.vegetation.soilvar.Owr) / (ground.STATVAR.vegetation.soilvar.Ows - ground.STATVAR.vegetation.soilvar.Owr))^(ground.STATVAR.vegetation.soilvar.n /(1-ground.STATVAR.vegetation.soilvar.n)-1))^(1/ground.STATVAR.vegetation.soilvar.n); %m
                
                % Actual water content is some fraction of saturation. These are only used for soil
                % thermal properties and phase change. Note the inconsistency with the use of soil
                % water in the bucket model to calculate the soil wetness factor.
                
                satfrac = 0.2;
                if (ground.STATVAR.vegetation.soilvar.tsoi(i) > 273.15)
                    ground.STATVAR.vegetation.soilvar.h2osoi_ice(i) = 0;
                    %         ground.STATVAR.vegetation.soilvar.h2osoi_liq(i) = satfrac * h2osoi_sat;
                else
                    %         ground.STATVAR.vegetation.soilvar.h2osoi_liq(i) = 0;
                    ground.STATVAR.vegetation.soilvar.h2osoi_ice(i) = satfrac * h2osoi_sat;
                end
                
            end  
            
        end
    end
end