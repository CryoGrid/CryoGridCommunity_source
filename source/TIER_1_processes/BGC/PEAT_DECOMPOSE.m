%========================================================================
% CryoGrid TIER1 library class PEAT_DECOMPOSE, containing functions related to peat decompostion in BGC_Frolking_peat
% S. Westermann, November 2021
%========================================================================

classdef PEAT_DECOMPOSE < BASE
    
    
    methods
        
        function peat = temp_modifier(peat)
            %peat.STATVAR.tempModifier = double(peat.STATVAR.T>= -10.0) .* exp(308.56 .* (1.0/40.02 - 1.0./(peat.STATVAR.T + 40.02)));

           %peat.STATVAR.tempModifier = peat.STATVAR.tempModifier .* double(peat.STATVAR.T >= 0);
           
           peat.STATVAR.tempModifier = double(peat.STATVAR.T >= 0) .* peat.PARA.Q10.^(peat.STATVAR.T ./ 10);
        end
        
        function peat = water_modifier(peat)
            peat.STATVAR.waterModifier = peat.STATVAR.vol_water .*0;
            range = peat.STATVAR.vol_water > peat.PARA.fieldCapacity;
            peat.STATVAR.waterModifier(range,1) = 1.0-(1.0-0.025).* ((peat.STATVAR.vol_water(range,1) - peat.PARA.fieldCapacity)./(1.0-peat.PARA.fieldCapacity)).^3.0;
            range = peat.STATVAR.vol_water >= 0.00 & peat.STATVAR.vol_water <= peat.PARA.fieldCapacity; %CHECK THIS!, should be 0.01 accoriding to Chaudhary 2017?
            peat.STATVAR.waterModifier(range,1) = 1.0-((peat.PARA.fieldCapacity - peat.STATVAR.vol_water(range,1))./peat.PARA.fieldCapacity).^5.0;%4.88);//IWRO 5 to 2//4.82
            %peat.STATVAR.tempModifier = 0.1;
        end
        
        function peat = peat_decompose(peat)
            peat.STATVAR.cato = peat.PARA.initialDecomposition.*(peat.STATVAR.total_peat./peat.STATVAR.totalpeatC_originalMass).^peat.PARA.decompo ;  %0.05
            peat.STATVAR.cato(isnan(peat.STATVAR.cato)) = 0; %if totalpeatC_originalMass = 0 
            peat.STATVAR.catm  = peat.STATVAR.cato .* peat.STATVAR.tempModifier .* peat.STATVAR.waterModifier' .* peat.PARA.decompose_timestep; 
             
            peat.TEMP.d_layerThick = peat.TEMP.d_layerThick - peat.STATVAR.catm.*peat.STATVAR.total_peat ./ peat.PARA.bulkDensity; % .* peat.CONST.mtocm; %add changes in bulk density
            peat.TEMP.d_organic = peat.TEMP.d_organic - peat.STATVAR.catm .* peat.STATVAR.total_peat ./ peat.CONST.organicDensity;
            
            
            peat.STATVAR.total_peat = peat.STATVAR.total_peat - peat.STATVAR.catm.*peat.STATVAR.total_peat;
            peat.STATVAR.layerThick = (peat.STATVAR.total_peat./peat.PARA.bulkDensity); %Original must be wrong??? NPP in kg/m2, buld density in kg/m3  .* peat.CONST.mtocm; % g/cm2 ./ cm3/g, then convert with m to cm

        end
        
        function peat = peat_decompose_Frolking(peat)
            
            %if ~isempty(peat.STATVAR.total_peat_PFT)
            peat.STATVAR.cato = repmat(peat.PARA.initial_decomposability, size(peat.STATVAR.total_peat_PFT,1),1).*(peat.STATVAR.total_peat_PFT./peat.STATVAR.totalpeatC_originalMass);  %0.05
%             else
%                 peat.STATVAR.cato=[];
%             end
            peat.STATVAR.cato(isnan(peat.STATVAR.cato)) = 0; %if totalpeatC_originalMass = 0 
            
            %ERC MODEL; REMOVE LATER
            peat.STATVAR.cato = peat.STATVAR.cato ./10;
            %END REMOVE
            
            peat.STATVAR.catm  = peat.STATVAR.cato .* repmat(peat.STATVAR.tempModifier,1,size(peat.PARA.initial_decomposability,2)) .* ...
                repmat(peat.STATVAR.waterModifier,1,size(peat.PARA.initial_decomposability,2)) .* peat.PARA.BGC_timestep; 
             
            peat.TEMP.d_layerThick = peat.TEMP.d_layerThick - sum(peat.STATVAR.catm.*peat.STATVAR.total_peat_PFT,2) ./ peat.PARA.bulkDensity; % .* peat.CONST.mtocm; %add changes in bulk density
            peat.TEMP.d_organic = peat.TEMP.d_organic - sum(peat.STATVAR.catm .* peat.STATVAR.total_peat_PFT, 2) ./ peat.CONST.organicDensity;
            
            
            peat.STATVAR.total_peat_PFT = peat.STATVAR.total_peat_PFT - peat.STATVAR.catm.*peat.STATVAR.total_peat_PFT;
            peat.STATVAR.total_peat = sum(peat.STATVAR.total_peat_PFT,2);
            peat.STATVAR.layerThick = (peat.STATVAR.total_peat./peat.PARA.bulkDensity); %Original must be wrong??? NPP in kg/m2, buld density in kg/m3  .* peat.CONST.mtocm; % g/cm2 ./ cm3/g, then convert with m to cm

        end
        

    end
end