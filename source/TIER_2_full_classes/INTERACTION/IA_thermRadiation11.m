%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction between two GROUND
% classes without water cycle
% contains function for SNOW CHILD phase 
% S. Westermann, October 2020
%========================================================================

classdef IA_thermRadiation11 <  IA_BASE 
    
    methods
        
        function get_boundary_condition_m(ia_heat, tile)
            stratigraphy1 = ia_heat.PREVIOUS;
            stratigraphy2 = ia_heat.NEXT;
            Isource1 = stratigraphy1.PARA.epsilon .* stratigraphy1.CONST.sigma.*(stratigraphy1.STATVAR.T(end) + stratigraphy1.CONST.Tmfw).^4;
            Isource2 = stratigraphy2.PARA.epsilon .* stratigraphy1.CONST.sigma.*(stratigraphy2.STATVAR.T(1) + stratigraphy1.CONST.Tmfw).^4;
            %residual = (stratigraphy2.PARA.epsilon -1) .* LW_out_strat1 + (stratigraphy1.PARA.epsilon -1) .* LW_out_strat2;
            LW_out_strat1 =  - Isource1 + Isource1 .* (1 - stratigraphy2.PARA.epsilon) .* stratigraphy1.PARA.epsilon ./ (1 - (1 - stratigraphy1.PARA.epsilon) .* (1 - stratigraphy2.PARA.epsilon));
            LW_out_strat2 =  - Isource2 + Isource2 .* (1 - stratigraphy1.PARA.epsilon) .* stratigraphy2.PARA.epsilon ./ (1 - (1 - stratigraphy1.PARA.epsilon) .* (1 - stratigraphy2.PARA.epsilon));
            
            stratigraphy1.TEMP.d_energy(end) = stratigraphy1.TEMP.d_energy(end) + (LW_out_strat1 - LW_out_strat2) .* stratigraphy1.STATVAR.area(end);
            stratigraphy2.TEMP.d_energy(1) = stratigraphy2.TEMP.d_energy(1) + (LW_out_strat2 - LW_out_strat1) .* stratigraphy2.STATVAR.area(1);
            
            stratigraphy1.TEMP.d_energy(end) = stratigraphy1.TEMP.d_energy(end) + stratigraphy1.PARA.exchangeAir .* tile.FORCING.TEMP.wind .* (-stratigraphy1.STATVAR.T(end) + tile.FORCING.TEMP.Tair);
            stratigraphy2.TEMP.d_energy(1) = stratigraphy2.TEMP.d_energy(1) + stratigraphy1.PARA.exchangeAir .* tile.FORCING.TEMP.wind .* (-stratigraphy2.STATVAR.T(1) + tile.FORCING.TEMP.Tair);
        end
         
%         %SNOW
%         function get_IA_CHILD_boundary_condition_u(ia_heat, tile)
%             get_boundary_condition_HEAT_IA_CHILD(ia_heat);
%         end
           
    end
end