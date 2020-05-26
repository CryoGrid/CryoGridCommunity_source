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
            
            
%           Flux aus canopy in Boden
            stratigraphy1.TEMP.F_lb = -stratigraphy1.STATVAR.vegetation.mlcanopyinst.gsoi;
            stratigraphy2.TEMP.F_ub = -stratigraphy1.TEMP.F_lb;
            
            %% function get_boundary_condition_water_m(ia_heat_vegetation)
            water2ground = stratigraphy1.STATVAR.vegetation.mlcanopyinst.qflx_prec_grnd_rain; % mm H20 / s
            stratigraphy2.TEMP.F_ub_water = water2ground / 1000.; % m H20 / s
            
            
            stratigraphy2.TEMP.dwc_dt = stratigraphy2.STATVAR.T .* 0;
            stratigraphy2.TEMP.dwc = stratigraphy2.STATVAR.T .* 0;
            
            stratigraphy1.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* stratigraphy1.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
            stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
            
            % wird in get_derivatives_prognostic(ground) dann ground übergeben, hier in fluxes geschrieben um der ordnung zu folgen und negativ
            stratigraphy2.TEMP.F_ub_water = -stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi-stratigraphy1.STATVAR.vegetation.soilvar.transp_per_layer(1) + stratigraphy2.TEMP.F_ub_water;
            stratigraphy2.TEMP.F_m_water(2:3) = -stratigraphy1.STATVAR.vegetation.soilvar.transp_per_layer(2)/2;
            stratigraphy2.TEMP.F_m_water(4:9) = -stratigraphy1.STATVAR.vegetation.soilvar.transp_per_layer(3)/6;
            stratigraphy2.TEMP.F_m_water(10:end) = 0;
        end
        
        function finalize_boundary_condition_water_m(ia_heat_vegetation,grid) % --> sollte in compute_diagnostics, ganz am Schluss ausgeführt werden im neuen Ground module. 
            
            stratigraphy1 = ia_heat_vegetation.PREVIOUS;
            stratigraphy2 = ia_heat_vegetation.NEXT;

            stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_vol(1) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.water,stratigraphy1.STATVAR.vegetation.soilvar.zi(1),'nearest');
            stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_vol(2) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.water,stratigraphy1.STATVAR.vegetation.soilvar.zi(2),'nearest');
            stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_vol(3) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.water,stratigraphy1.STATVAR.vegetation.soilvar.zi(3),'nearest');
            stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_ice(1) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.ice,stratigraphy1.STATVAR.vegetation.soilvar.zi(1),'nearest'); 
            stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_ice(2) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.ice,stratigraphy1.STATVAR.vegetation.soilvar.zi(2),'nearest');
            stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_ice(3) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.ice,stratigraphy1.STATVAR.vegetation.soilvar.zi(3),'nearest');
            
%             stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_vol(1) = stratigraphy2.STATVAR.water(1);
%             stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_vol(2) = (stratigraphy2.STATVAR.water(2)+stratigraphy2.STATVAR.water(3))/2;
%             stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_vol(3) = (stratigraphy2.STATVAR.water(4)+stratigraphy2.STATVAR.water(5)+stratigraphy2.STATVAR.water(6)+stratigraphy2.STATVAR.water(7)+stratigraphy2.STATVAR.water(8)+stratigraphy2.STATVAR.water(9))/6; 
%             stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_ice(1) = stratigraphy2.STATVAR.ice(1);
%             stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_ice(2) = (stratigraphy2.STATVAR.ice(2)+stratigraphy2.STATVAR.ice(3))/2;
%             stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_ice(3) = (stratigraphy2.STATVAR.ice(4)+stratigraphy2.STATVAR.ice(5)+stratigraphy2.STATVAR.ice(6)+stratigraphy2.STATVAR.ice(7)+stratigraphy2.STATVAR.ice(8)+stratigraphy2.STATVAR.ice(9))/6; 
                  
            % Simone: Set ground.STATVAR.vegetation.soilvar.t_soisno to actual CG soil temperature            
            stratigraphy1.STATVAR.vegetation.soilvar.t_soisno(1) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.T+273.15,stratigraphy1.STATVAR.vegetation.soilvar.zi(1),'nearest');
            stratigraphy1.STATVAR.vegetation.soilvar.t_soisno(2) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.T+273.15,stratigraphy1.STATVAR.vegetation.soilvar.zi(2),'nearest');
            stratigraphy1.STATVAR.vegetation.soilvar.t_soisno(3) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.T+273.15,stratigraphy1.STATVAR.vegetation.soilvar.zi(3),'nearest');
    
        end
    end
end