classdef PEAT_ACCUMULATION < BASE
    
    
    methods
        
%         %called only once in finalize_init
%         function peat = readingNppData(peat)  
%             data = load('anpp_whyme.mat'); %loading ANPP data (KgC/m-2/year)
%             ANPP = data.data(:,30);
%             anpp = ANPP(1:peat.PARA.noy,1);%data
%             
%             anppC = repmat(anpp,[1 peat.PARA.noi]); %creating repetitive npp values based on noi
%             app_dim1 = size(anppC,1); % first dimension of anpp matrix
%             app_dim2 = size(anppC,2); % second dimension of anpp matrix
%             anpp_dim = app_dim1 * app_dim2; % multiplying both dimensions
%             anpp = reshape (anppC,anpp_dim,1);
%             
%             boundry_condition = 200;%%23.5; %%25;% added 60 kg C for the boundary condition
%             Anpp = anpp;%[boundry_condition;anpp]; % appendding boundry conditions before actual anpp data
%             peat.STATVAR.Anpp =Anpp(1:peat.PARA.tyears);
%         end
%         
%         %called only once in finalize_init
%         function peat = distributing_NPP(peat)
%             % %foliar projective cover
%             ss = 110;%109 + rept;%110
%             rand('seed', ss); %#ok<RAND>
%             
%             peat.STATVAR.fpcs = rand (1,peat.PARA.nop);
%             peat.STATVAR.fpcg = rand (1,peat.PARA.nop);
%             peat.STATVAR.fpcm = rand (1,peat.PARA.nop);
%             
%             peat.STATVAR.fpc = peat.STATVAR.fpcm+ peat.STATVAR.fpcs + peat.STATVAR.fpcg;
%             peat.STATVAR.fpcs = peat.STATVAR.fpcs./peat.STATVAR.fpc;
%             peat.STATVAR.fpcm = peat.STATVAR.fpcm./peat.STATVAR.fpc;
%             peat.STATVAR.fpcg = peat.STATVAR.fpcg./peat.STATVAR.fpc;
%             peat.STATVAR.fpc = peat.STATVAR.fpcm+ peat.STATVAR.fpcs + peat.STATVAR.fpcg;
%             peat.STATVAR.fpccm = peat.STATVAR.fpcm;
%             peat.STATVAR.fpccg = peat.STATVAR.fpcg;
%             peat.STATVAR.fpccs = peat.STATVAR.fpcs;
%             
%             % % Annual npp
%             peat.STATVAR.npp_shrubC = peat.PARA.npp_frac_shrub*(peat.STATVAR.Anpp*peat.STATVAR.fpcs); % npp of shrubs in every patch
%             peat.STATVAR.npp_mossC = peat.STATVAR.Anpp*peat.STATVAR.fpcm; % npp of mosses in every patch
%             peat.STATVAR.npp_graminoidC = peat.PARA.npp_frac_graminoid*(peat.STATVAR.Anpp*peat.STATVAR.fpcg);  % npp of grasses in every patch
%             peat.STATVAR.tnpp = peat.STATVAR.npp_shrubC+peat.STATVAR.npp_mossC+peat.STATVAR.npp_graminoidC;% total npp
%             
%             %%normalising npp
%             peat.STATVAR.npp_shrub_frac = (peat.STATVAR.npp_shrubC./peat.STATVAR.tnpp);
%             peat.STATVAR.npp_graminoid_frac = (peat.STATVAR.npp_graminoidC./peat.STATVAR.tnpp);
%             peat.STATVAR.npp_moss_frac= (peat.STATVAR.npp_mossC./peat.STATVAR.tnpp);
%             
%             for s=1:peat.PARA.nop
%                 peat.STATVAR.npp_shrub(:,s)  = peat.STATVAR.npp_shrub_frac(:,s).* peat.STATVAR.Anpp;
%                 peat.STATVAR.npp_graminoid(:,s)  = peat.STATVAR.npp_graminoid_frac(:,s).* peat.STATVAR.Anpp;
%                 peat.STATVAR.npp_moss(:,s)  = peat.STATVAR.npp_moss_frac(:,s).* peat.STATVAR.Anpp;
%             end
%             
%             peat.STATVAR.tnpp1 = peat.STATVAR.npp_shrub+peat.STATVAR.npp_moss+peat.STATVAR.npp_graminoid;
%             
%             peat.STATVAR.npp_Shrub = peat.STATVAR.npp_shrub;
%             peat.STATVAR.npp_Moss = peat.STATVAR.npp_moss;
%             peat.STATVAR.npp_Graminoid = peat.STATVAR.npp_graminoid;
%             
%         end
        
%         function peat = get_peatC(peat)
%            
%             peat.STATVAR.peat_moss = peat.STATVAR.npp_moss(peat.TEMP.year_count+1,:);%  moss peat
%             peat.STATVAR.peat_shrub = peat.STATVAR.npp_shrub(peat.TEMP.year_count+1,:); % shrub peat
%             peat.STATVAR.peat_graminoid = peat.STATVAR.npp_graminoid(peat.TEMP.year_count+1,:); % graminoid peat
% 
%         end
        

        function peat = peat_accumulation_newCell(peat)
            
             %new_peat = (peat.STATVAR.peat_moss + peat.STATVAR.peat_shrub + peat.STATVAR.peat_graminoid) ./ peat.PARA.number_of_growing_days; % total peat, in g/cm2!!!
                         
             %assume constant NPP
             new_peat = peat.PARA.annual_NPP./ peat.PARA.number_of_growing_days .* peat.PARA.decompose_timestep;
             
             peat.STATVAR.total_peat = [new_peat; peat.STATVAR.total_peat];
             peat.STATVAR.totalpeatC_originalMass = [new_peat; peat.STATVAR.totalpeatC_originalMass]; %;  peat.STATVAR.totalpeatC_originalMass_old];
            
            peat.STATVAR.layerThick = (peat.STATVAR.total_peat./peat.PARA.bulkDensity);
            
            peat.TEMP.d_layerThick = [new_peat ./ peat.PARA.bulkDensity; peat.TEMP.d_layerThick];
            peat.TEMP.d_organic = [ new_peat ./ peat.CONST.organicDensity ; peat.TEMP.d_organic];
            
        end
        
        function peat = peat_accumulation(peat)
            
             %new_peat = (peat.STATVAR.peat_moss + peat.STATVAR.peat_shrub + peat.STATVAR.peat_graminoid)  ./ peat.PARA.number_of_growing_days; % total peat, in kg/m2 ??
          
             %assume constant NPP
             new_peat = peat.PARA.annual_NPP./ peat.PARA.number_of_growing_days .* peat.PARA.decompose_timestep;
             
             
             peat.STATVAR.total_peat(1,1) = peat.STATVAR.total_peat(1,1) + new_peat;
             peat.STATVAR.totalpeatC_originalMass(1,1) = peat.STATVAR.totalpeatC_originalMass(1,1) + new_peat; %;  peat.STATVAR.totalpeatC_originalMass_old];
             
             peat.STATVAR.layerThick = (peat.STATVAR.total_peat./peat.PARA.bulkDensity); %CONVERSION MUST BE WRONG, SEE DECOMPOSITION.* peat.CONST.mtocm; % g/cm2 ./ cm3/g, then convert with m to cm
             
             peat.TEMP.d_layerThick(1,1) = peat.TEMP.d_layerThick(1,1) + new_peat ./ peat.PARA.bulkDensity;
             peat.TEMP.d_organic(1,1) = peat.TEMP.d_organic(1,1) + new_peat ./ peat.CONST.organicDensity;
             
        end
        
        function peat = updatebulkD(peat)
            
            peat.STATVAR.bulkDensity = (peat.STATVAR.organic> 0).*peat.PARA.minbulkDensity+(peat.PARA.diffbulkDensity .*(1/(1+exp(-(40*(1-peat.STATVAR.massRemain(1))-34))))) + (peat.STATVAR.organic <= 0).* peat.PARA.mineral_bulkDensity;
            
        end
        
        
        %not sure these functions are required at some stage
        
        %         function peat = random_topo (peat)
%             reptt = 777;
%             %%generating random topography
%             rand('seed', 110+reptt); %#ok<RAND>
%             bcc = (peat.PARA.initial_topo-peat.PARA.initial_topo*peat.PARA.adjusting_SH); % making uneven peat
%             peat.STATVAR.total_peat (1,:) =  bcc + rand (1,peat.PARA.nop)*(peat.PARA.initial_topo-bcc);
%             peat.PARA.initial_topo = peat.STATVAR.total_peat (1,:);
%         end
%         function peat = get_boundary_condition_peat_u(peat)
%             
%             %use TIER1 functions, add NPP, etc
%             %if (peat.PARA.dayofyear ==1 && peat.PARA.year == 1 && peat.PARA.flag ==1)  %on the first day of the first year
%                 peat = readingNppData(peat);  %reading NPP data
%                 peat = distributing_NPP(peat); %distributing NPP data in different patches
%               %  peat.PARA.flag = 0;
% %             end
%             
%            % if(peat.PARA.dayofyear ==1 && peat.PARA.flag2 ==1)
%                 
%                 peat = get_peatC(peat);
% %                 if peat.PARA.year == 1, peat = random_topo (peat); end  %creating random topography
%                 %  e.g. set upper boundary T to 5 degree C
%             %end
%             
%             %Temp. and water in the upper boundary - daily timestep
% %             peat = get_ub_T_WaterAndSnow(peat);
%             
%             %initializing water and ice variable based on peat density and
%             %peat_accumulation
%           %  if (peat.PARA.dayofyear ==1 && peat.PARA.flag3 ==1)
% %                 peat = peat_accumulation(peat);
%            % end
%             
%         end
        
%         function peat = accumulate_peat(peat)
%             
%             peat.STATVAR.tpeat = peat.STATVAR.total_peat(1,:);
%             
%             peat.STATVAR.total_cumulative_peat = cumsum(peat.STATVAR.total_peat);
%             
%             if peat.PARA.dayofyear ==1
%                 peat.STATVAR.total_cumulative_peat(1,1:peat.PARA.nop)= peat.STATVAR.total_peat(1,1:peat.PARA.nop);%removing error of cumsum function on the first year
%             end
%             
%             peat.STATVAR.peat_depth = 0.0 + (peat.STATVAR.total_peat./peat.STATVAR.bulkDensity)* 100; % m to cm
%             peat.STATVAR.peatD_all = peat.STATVAR.peat_depth(peat.PARA.dayofyear,:)+peat.STATVAR.peatD_all;%total peat depth every year (cm)
%             
%             peat.STATVAR.cumulative_peat  = cumsum(peat.STATVAR.peat_depth);% total cumulative peat depth each year (cm)
%             
%             if peat.PARA.dayofyear == 1
%                 peat.STATVAR.cumulative_peat(1,1:peat.PARA.nop)= peat.STATVAR.peat_depth(1,1:peat.PARA.nop);
%             end
%             
%             peat.STATVAR.ccpeat = peat.STATVAR.cumulative_peat(end,:);
%             
%             peat.STATVAR.max_waterHoldingCapacity = peat.STATVAR.peat_depth* peat.PARA.porosity;%(cm)
%             
%             peat.STATVAR.mwater = cumsum(peat.STATVAR.max_waterHoldingCapacity);
%             
%             if peat.PARA.dayofyear == 1
%                 peat.STATVAR.mwater(1,(1:peat.PARA.nop))= peat.STATVAR.max_waterHoldingCapacity(1,(1:peat.PARA.nop));% fixing the last value of first row--error coming due to cumsum function
%             end
%             
%             peat.STATVAR.mx = peat.STATVAR.mwater(end,:)*10;% + peat.STATVAR.heightOfWater;  % total accumulative water holding capacity (mm)
%             
%         end
        
%         function peat = get_derivative_PEAT(peat)
%             peat.TEMP.C_derivative = peat.STATVAR.peat_moss(peat.PARA.year) + peat.STATVAR.peat_shrub(peat.PARA.year) + peat.STATVAR.peat_graminoid(peat.PARA.year); % total peat
%             peat.TEMP.C_derivative = peat.TEMP.C_derivative;%/365;
%             peat.PARA.t_old = peat.PARA.t;
%             %             peat.PARA.t;
%         end
%         
%         function peat = calculate_massRemain(peat)
%             
%             %             peat.STATVAR.massRemain(1) =
%             
%         end
        

      
    end
end