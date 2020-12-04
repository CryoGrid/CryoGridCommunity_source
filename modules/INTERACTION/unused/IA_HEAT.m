classdef IA_HEAT < matlab.mixin.Copyable 

     properties
        PREVIOUS
        NEXT
    end
    
    methods

        function  get_boundary_condition_m(ia_heat)
            stratigraphy1 = ia_heat.NEXT; %Vegetation            
            stratigraphy2 = ia_heat.PREVIOUS; %Snow 
            stratigraphy3 = ia_heat.NEXT.NEXT; %Ground
            
            flux = (stratigraphy2.STATVAR.T(end) - stratigraphy3.STATVAR.T(1)) .* stratigraphy2.STATVAR.thermCond(end) .* stratigraphy3.STATVAR.thermCond(1) ./...
                (stratigraphy2.STATVAR.thermCond(end).* stratigraphy3.STATVAR.layerThick(1)./2 + stratigraphy3.STATVAR.thermCond(1).* stratigraphy2.STATVAR.layerThick(end)./2 );
  
            stratigraphy1.TEMP.F_lb = -stratigraphy1.STATVAR.vegetation.mlcanopyinst.gsoi;
            stratigraphy2.TEMP.F_ub = -stratigraphy1.TEMP.F_lb; 
            
            % Added by Simone
            stratigraphy2.TEMP.F_lb = -flux;
            stratigraphy3.TEMP.F_ub = flux;
  
            % Water fluxes between Vegetation and Ground
            water2ground = stratigraphy1.STATVAR.vegetation.mlcanopyinst.qflx_prec_grnd_rain; % mm H20 / s
            stratigraphy2.TEMP.F_ub_water = water2ground / 1000.; % m H20 / s
            
            stratigraphy1.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* stratigraphy1.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
            stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
            
            % M [kg/mol] * transp/etsoi [mol/m2s] / dichte H2O [kg/m3]
            % dwc_dt [m/s]
            
            % wird in get_derivatives_prognostic(ground) dann ground übergeben, hier in fluxes geschrieben um der ordnung zu folgen und negativ
            stratigraphy2.TEMP.F_ub_water = -stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi + stratigraphy2.TEMP.F_ub_water;
            
            stratigraphy3.TEMP.F_ub_water = -stratigraphy1.STATVAR.vegetation.soilvar.transp_per_layer(1) + stratigraphy2.TEMP.F_ub_water;
            stratigraphy3.TEMP.F_m_water(2:3) = -stratigraphy1.STATVAR.vegetation.soilvar.transp_per_layer(2)/2;
            stratigraphy3.TEMP.F_m_water(4:9) = -stratigraphy1.STATVAR.vegetation.soilvar.transp_per_layer(3)/6;
            stratigraphy3.TEMP.F_m_water(10:end) = 0;
            
        end
        
        function finalize_boundary_condition_water_m(ia_heat,grid) % --> sollte in compute_diagnostics, ganz am Schluss ausgeführt werden im neuen Ground module. 
        end
        
    end
end
