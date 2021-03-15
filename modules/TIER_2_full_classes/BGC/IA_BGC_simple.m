%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction between two GROUND
% classes without water cycle
% contains function for SNOW CHILD phase 
% S. Westermann, October 2020
%========================================================================

classdef IA_BGC_simple <  IA_BGC
    
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
            depths_ground =  cumsum([0; ia_BGC.GROUND.STATVAR.layerThick]);
            depths_BGC =  cumsum([0; ia_BGC.BGC.STATVAR.layerThick]);
            %set everything to zero 
            ia_BGC.BGC.STATVAR.vol_water = ia_BGC.BGC.STATVAR.layerThick .*0;
            ia_BGC.BGC.STATVAR.vol_mineral = ia_BGC.BGC.STATVAR.layerThick .*0;
            ia_BGC.BGC.STATVAR.porosity = ia_BGC.BGC.STATVAR.layerThick .*0;
            ia_BGC.BGC.STATVAR.T = ia_BGC.BGC.STATVAR.layerThick .*0;
            ia_BGC.BGC.STATVAR.field_capacity = ia_BGC.BGC.STATVAR.layerThick .*0;
            vol_water_ground = ia_BGC.GROUND.STATVAR.water ./ ia_BGC.GROUND.STATVAR.layerThick ./ ia_BGC.GROUND.STATVAR.area;
            vol_mineral_ground = ia_BGC.GROUND.STATVAR.mineral ./ ia_BGC.GROUND.STATVAR.layerThick ./ ia_BGC.GROUND.STATVAR.area;    
            porosity_ground = 1- (ia_BGC.GROUND.STATVAR.mineral + ia_BGC.GROUND.STATVAR.organic) ./ ia_BGC.GROUND.STATVAR.layerThick ./ ia_BGC.GROUND.STATVAR.area;    
            
            
            overlap = get_overlap_cells(ia_BGC, depths_BGC, depths_ground);
            ia_BGC.BGC.STATVAR.overlap = overlap;
            
            for i=1:size(overlap,1)
                ia_BGC.BGC.STATVAR.T(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.T(overlap(i, 1),1) + overlap(i,3) .* ia_BGC.GROUND.STATVAR.T(overlap(i, 2),1); 
                ia_BGC.BGC.STATVAR.vol_water(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.vol_water(overlap(i, 1),1) + overlap(i,3) .* vol_water_ground(overlap(i, 2),1); 
                ia_BGC.BGC.STATVAR.vol_mineral(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.vol_mineral(overlap(i, 1),1) + overlap(i,3) .* vol_mineral_ground(overlap(i, 2),1); 
                ia_BGC.BGC.STATVAR.porosity(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.porosity(overlap(i, 1),1) + overlap(i,3) .* porosity_ground(overlap(i, 2),1); 
                ia_BGC.BGC.STATVAR.field_capacity(overlap(i, 1),1) = ia_BGC.BGC.STATVAR.field_capacity(overlap(i, 1),1) + overlap(i,3) .* ia_BGC.GROUND.STATVAR.field_capacity(overlap(i, 2),1);
            end
            
        end
           
    end
end