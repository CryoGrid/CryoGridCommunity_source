classdef PEAT_DECOMPOSE < BASE
    
    properties
        
    end
    
    methods
        
          %COMMENT Sebatsian: is this really only done for the first cell??
        function peat = temp_modifier(peat)
            
%             for i = 1:length(peat.STATVAR.T)
%             if (peat.STATVAR.T(i,1) >= -10.0)
%                 
%                 peat.STATVAR.tempModifier(i,1) = exp(308.56 .* (1.0/40.02 - 1.0/(peat.STATVAR.T(i,1)+40.02))); %% gtemp
%             else
%                 peat.STATVAR.tempModifier(i,1) = 0.0;
%             end
%             end
            peat.STATVAR.tempModifier = double(peat.STATVAR.T>= -10.0) .* exp(308.56 .* (1.0/40.02 - 1.0./(peat.STATVAR.T + 40.02)));
        end
        
        function peat = water_modifier(peat)
            
%             for i = 1:peat.PARA.year
%             %added Sebastian, waterContent is always 0 now:
%             peat.STATVAR.waterContent(i,1) = peat.STATVAR.water (i,1) ./ (peat.STATVAR.peat_depth(i,1));%.* peat.STATVAR.area (i,1)- peat.STATVAR.XwaterIce (i,1));
%             %this is the volumetric water content
%             
%             if (peat.STATVAR.waterContent(i,1) > peat.PARA.fieldCapacity)  % soil water modifier
%                 
%                 peat.STATVAR.waterModifier(i,1) = 1.0-(1.0-0.025).* ((peat.STATVAR.waterContent(i,1) -peat.PARA.fieldCapacity)/(1.0-peat.PARA.fieldCapacity)).^5.0;
%                 
%             elseif (peat.STATVAR.waterContent(i,1) < 0.1)% && wtp_patch[0][p2] < -400)
%                 
%                 
%                 peat.STATVAR.waterModifier(i,1)  = 0.6;
%             else
%                 
%                 peat.STATVAR.waterModifier(i,1)  = 1.0-((peat.PARA.fieldCapacity-peat.STATVAR.waterContent(i,1))/peat.PARA.fieldCapacity).^5.0;%4.88);//IWRO 5 to 2//4.82
%             end
%             end
            
            range = peat.STATVAR.vol_water > peat.STATVAR.field_capacity;
            peat.STATVAR.waterModifier(range) = 1.0-(1.0-0.025).* ((peat.STATVAR.vol_water(range) - peat.STATVAR.field_capacity(range))./(1.0-peat.STATVAR.field_capacity(range))).^5.0;
            range = peat.STATVAR.vol_water < 0.01;
            peat.STATVAR.waterModifier(range) = 0.6;
            range = peat.STATVAR.vol_water >= 0.01 & peat.STATVAR.vol_water <= peat.STATVAR.field_capacity; %CHECK THIS!, should be 0.01 accoriding to Chaudhary 2017?
            peat.STATVAR.waterModifier(range) = max(0.6, 1.0-((peat.STATVAR.field_capacity(range) - peat.STATVAR.vol_water(range))./peat.STATVAR.field_capacity(range)).^5.0);%4.88);//IWRO 5 to 2//4.82
            
        end
        
        function peat = peat_decompose(peat)
            peat.STATVAR.cato = peat.PARA.initialDecomposition*(peat.STATVAR.total_peat./peat.STATVAR.totalpeatC_originalMass).^peat.PARA.decompo ;  %0.05
            peat.STATVAR.cato(isnan(peat.STATVAR.cato)) = 0; %if totalpeatC_originalMass = 0 
            peat.STATVAR.catm  = peat.STATVAR.cato .* peat.STATVAR.tempModifier .* peat.STATVAR.waterModifier; 
            peat.STATVAR.total_peat = peat.STATVAR.total_peat - peat.STATVAR.catm.*peat.STATVAR.total_peat;
            
          
        end
        
  
    end
end