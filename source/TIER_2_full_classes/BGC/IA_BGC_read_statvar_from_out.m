%========================================================================
% CryoGrid INTERACTION (IA) for GROUND_..._BGC classes with
% GROUND_store_flip_flop classes
% S. Westermann, November 2021
%========================================================================

classdef IA_BGC_read_statvar_from_out <  IA_BGC
    
    properties
        GROUND
        BGC
    end
    
    methods
        
        function send_BGC_variables(ia_BGC, tile)
            %do nothing at this point - for fully coupled runs, this must
            %must include transfer of mass to and within
            %the physics stratigraphy
        end
        
        function get_ground_variables(ia_BGC, tile)
            %map the physical variables required by the BGC module to the BGC grid
            %depths_ground =  cumsum([0; ia_BGC.GROUND.STATVAR.layerThick]);
            depths_ground =  cumsum([0; ia_BGC.GROUND.STATVAR.layerThick - ia_BGC.GROUND.STATVAR.XwaterIce ./ ia_BGC.GROUND.STATVAR.area]);
            depths_BGC =  cumsum([0; ia_BGC.BGC.STATVAR.layerThick]);
            %set everything to zero 
            ia_BGC.BGC.STATVAR.vol_water = ia_BGC.BGC.STATVAR.layerThick .*0;
%             ia_BGC.BGC.STATVAR.vol_mineral = ia_BGC.BGC.STATVAR.layerThick .*0;
%             ia_BGC.BGC.STATVAR.porosity = ia_BGC.BGC.STATVAR.layerThick .*0;
            ia_BGC.BGC.STATVAR.T = ia_BGC.BGC.STATVAR.layerThick .*0;
            ia_BGC.BGC.STATVAR.field_capacity = ia_BGC.BGC.STATVAR.layerThick .*0;
           % vol_water_ground = ia_BGC.GROUND.STATVAR.water ./ ia_BGC.GROUND.STATVAR.layerThick ./ ia_BGC.GROUND.STATVAR.area;
            vol_water_ground = ia_BGC.GROUND.STATVAR.waterIce ./ (ia_BGC.GROUND.STATVAR.layerThick .* ia_BGC.GROUND.STATVAR.area - ia_BGC.GROUND.STATVAR.XwaterIce - ia_BGC.GROUND.STATVAR.organic - ia_BGC.GROUND.STATVAR.mineral);
%             vol_mineral_ground = ia_BGC.GROUND.STATVAR.mineral ./ ia_BGC.GROUND.STATVAR.layerThick ./ ia_BGC.GROUND.STATVAR.area;    
%             porosity_ground = 1- (ia_BGC.GROUND.STATVAR.mineral - ia_BGC.GROUND.STATVAR.organic) ./ ia_BGC.GROUND.STATVAR.layerThick ./ ia_BGC.GROUND.STATVAR.area;    
            
            
            overlap = get_overlap_cells(ia_BGC, depths_BGC, depths_ground);
            ia_BGC.BGC.STATVAR.overlap = overlap;
            
            for i=1:size(overlap,1)
                ia_BGC.BGC.STATVAR.T(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.T(overlap(i, 1),1) + overlap(i,3) .* ia_BGC.GROUND.STATVAR.T(overlap(i, 2),1); 
                ia_BGC.BGC.STATVAR.vol_water(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.vol_water(overlap(i, 1),1) + overlap(i,3) .* vol_water_ground(overlap(i, 2),1); 
                %a_BGC.BGC.STATVAR.vol_mineral(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.vol_mineral(overlap(i, 1),1) + overlap(i,3) .* vol_mineral_ground(overlap(i, 2),1); 
                %ia_BGC.BGC.STATVAR.porosity(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.porosity(overlap(i, 1),1) + overlap(i,3) .* porosity_ground(overlap(i, 2),1); 
                ia_BGC.BGC.STATVAR.field_capacity(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.field_capacity(overlap(i, 1),1) + overlap(i,3) .* ia_BGC.GROUND.STATVAR.field_capacity(overlap(i, 2),1);
            end
            
        end
        
        function peat_depth = get_peat_depth(ia_BGC, tile)
            %INCLUDE XICE, USE OVERSATURATION IN GROUND CELL AND COMPUTE
            %OVERLAP
            peat_depth = sum(ia_BGC.BGC.STATVAR.layerThick,1);
        end
           
        function get_NPP_Frolking(ia_BGC, tile)
            T = ia_BGC.GROUND.STATVAR.T(1,1);
            temp_modifier = max(0,min(1,min((T - ia_BGC.BGC.PARA.start_PhotSyn_T)./(ia_BGC.BGC.PARA.start_fullPhotSyn_T - ia_BGC.BGC.PARA.start_PhotSyn_T), ...
                (T - ia_BGC.BGC.PARA.end_PhotSyn_T)./ (ia_BGC.BGC.PARA.end_fullPhotSyn_T - ia_BGC.BGC.PARA.end_PhotSyn_T))));
            %disp([T temp_modifier])
            radiation_modifier = ia_BGC.GROUND.STATVAR.Sin;
            ET_modifier = max(0, ia_BGC.GROUND.STATVAR.Qe ./ ia_BGC.GROUND.STATVAR.Qe_pot);
            
            ia_BGC.BGC.TEMP.GPP = max(0, temp_modifier .* radiation_modifier .* ET_modifier);
            
            %get water table depth
            saturation = ia_BGC.GROUND.STATVAR.waterIce ./ (ia_BGC.GROUND.STATVAR.layerThick .* ia_BGC.GROUND.STATVAR.area - ia_BGC.GROUND.STATVAR.XwaterIce - ia_BGC.GROUND.STATVAR.mineral - ia_BGC.GROUND.STATVAR.organic);
            ia_BGC.BGC.TEMP.water_table_depth = sum(ia_BGC.GROUND.STATVAR.layerThick(1:find(saturation>0.95,1)-1,1),1) - ia_BGC.GROUND.STATVAR.XwaterIce(1,1)./ ia_BGC.GROUND.STATVAR.area(1,1);
            ia_BGC.BGC.TEMP.water_table_depth = min(ia_BGC.BGC.TEMP.water_table_depth, 1.5);
        end
        
        
        function add_grid_cell(ia_BGC, tile)
            ia_BGC.BGC.STATVAR.BGC_overlap_vector = [1; ia_BGC.BGC.STATVAR.BGC_overlap_vector];
            
        end
        
        function regrid_stratigraphy(ia_BGC, tile)
            
        end
    end
end