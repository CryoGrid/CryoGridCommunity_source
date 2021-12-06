%========================================================================
% CryoGrid INTERACTION (IA) class for biogeochemistry BGC classes with
% physical GROUND classes
% S. Westermann, November 2021
%========================================================================

classdef IA_BGC <  IA_BASE
    
    properties

    end
    
    methods
        
        %service functions for BGC IA classes
        
        function overlap = get_overlap_cells(ia_BGC, cell_1, cell_2)
            
            overlap = [];
            
            if size(cell_1,1) > 1 && size(cell_2,1) > 1
                for i1=1:size(cell_1,1)-1
                    i2=1;
                    a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                    
                    while a <= 0 && i2 < size(cell_2,1)-1
                        %overlap2(i1,i2) = a;
                        i2 = i2+1;
                        a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                    end
                    if a>0
                        overlap = [overlap;  [i1  i2 a ./ ia_BGC.BGC.STATVAR.layerThick(i1,1) a ./ ia_BGC.GROUND.STATVAR.layerThick(i2,1)]];
                    end
                    
                    i2_start = i2;
                    while a > 0 && i2 < size(cell_2,1)-1
                        
                        i2 = i2+1;
                        a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                        if a>0
                            overlap = [overlap;  [i1  i2 a ./ ia_BGC.BGC.STATVAR.layerThick(i1,1) a ./ ia_BGC.GROUND.STATVAR.layerThick(i2,1)]];
                        end
                    end
                    i2 = i2_start;
                end
            end
        end

           
    end
end