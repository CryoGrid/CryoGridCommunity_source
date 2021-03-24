%========================================================================
% CryoGrid TIER1 library class for functions related to the CLM4 biogeochemistry

% S. Westermann, October 2020
%========================================================================

classdef BGC_CLM4 < BASE
    
    methods
        function [ GRID, Crate, gasRate, totalCO2rate] = aerobicRate( GRID, PARA, T, Crate, gasRate)
            %AERODECOMP Simulates the decomposition of soil organic matter and
            %associated consumption and production of gases.
            
            R = 8.314;   %[J / mol K] Universal gas constant
            
            % Read conditions:
            
            %Find the soil cells in the BGC domain:
            decompDomain = GRID.soil.cT_domain(GRID.BGC.cT_domain);
            
            %Calculate the number of cells that contain snow:
            nSoilCells = sum(double(decompDomain));
            T = T(GRID.soil.cT_domain & GRID.BGC.cT_domain);  %Temperature in the soil cells in the BGC domain.
            
            %Volume fractions:
            thetaSat = 1-GRID.soil.cT_organic(GRID.BGC.cT_domain(GRID.soil.cT_domain))-GRID.soil.cT_mineral(GRID.BGC.cT_domain(GRID.soil.cT_domain));                 %[-]
            
            %Volume of air and water in each grid cell [m^3/m2 = m]
            totVolume = thetaSat.*GRID.general.K_delta(GRID.BGC.cT_domain & GRID.soil.cT_domain);                %[m]
            
            %Read the total masses of O2 and CO2:
            m_gases = GRID.BGC.m_gases(decompDomain,1:2);
            
            %Calculate the overall concentrations in the pores:
            O2PerPorV = m_gases(:,1)./totVolume;  %[g/m^3]
            
            Q10 = zeros(size(GRID.BGC.C_mass));
            %Calculate the temperature effect:
            for i =1:size(Q10,2)
                Q10(:,i) = PARA.BGC.Q10pools(i).^( (T-PARA.BGC.T_ref(i))/10 );       %[-]
            end
            
            % Q10(T<0,:)= 0;
            
            %%
            %Calculate the O2-limitation on the rate:
            kappaO2 = O2PerPorV./(PARA.BGC.HT_halfRO2 + O2PerPorV);       %[-]
            
            %% Scale the rate down in proportion to the ice-filled porous volume:
            kappaIce = 1- GRID.BGC.thetaIceSoil./GRID.BGC.thetaSatSoil;
            %% Calculate the CO2-limitation on the rate:
            %Calculate the cut-off solubility-weighted concentration of CO2:
            % cCO2cutOff = PARA.BGC.alpha_CO2 * 5.36*10^5 ./(T+273.15);          %[g/m^3]
            % cCO2 = GRID.BGC.c(decompDomain,2);
            %
            % kappa_CO2 = zeros(nSoilCells,1);
            % belowLimit = cCO2< cCO2cutOff;
            % if any(belowLimit)
            %     kappa_CO2(belowLimit) = (1-(cCO2(belowLimit)./cCO2cutOff(belowLimit)).^4);     %[-]
            % end
            kappa_CO2 = 1;
            
            %%
            %Calculate the water-potential limitation on the rate:
            Psi = GRID.BGC.SoilwaterPot;
            PsiMin = PARA.BGC.PsiMin;
            PsiMax = PARA.BGC.PsiMax;
            
            %TODO: Implement reduction in water potenital during drought.
            kappa_water = zeros(nSoilCells,1);                          %[-]
            kappa_water(Psi > PsiMin & Psi< PsiMax) = ...
                log10(PsiMin./Psi(Psi > PsiMin & Psi< PsiMax))/log10(PsiMin/PsiMax);%[-]
            kappa_water(Psi >= PARA.BGC.PsiMax) = 1;                    %[-]
            % Store water rate limiter for later:
            GRID.BGC.kappa_water = kappa_water;
            
            % Calculate the threshold pressure.
            P_threshold = GRID.BGC.waterP(decompDomain)*PARA.BGC.modelledGasFraction;
            
            %Calculate the partial pressures of all modelled gases:
            sumPartials = R*(T+273.15).*(GRID.BGC.c(decompDomain,:)*[32; 44; 16].^(-1));
            
            kappaPressure = max(1-(sumPartials./P_threshold).^4, 0);
            
            %Calculate non-pool-specific rate coefficients:
            KAPPA_nps = kappaO2.* kappa_CO2 .* kappa_water.*kappaIce.*kappaPressure;
            % %Calculate non-pool-specific rate coefficients, ignoring CO2 deprecation:
            % KAPPA_nps = kappa_O2 .* kappa_water;
            
            %Calculate pool-specific rates:
            poolSpecRates = Q10.*(ones(nSoilCells,1)*PARA.BGC.k_pool);
            
            
            % Calculate the compound rate of C turnover
            compoundRates = KAPPA_nps*ones(1,size(GRID.BGC.C_mass,2)).*poolSpecRates;
            
            decompRates = GRID.BGC.C_mass .* compoundRates;
            
            % Calculate the translocation rates:
            translocOut = decompRates.*(1-ones(size(decompRates,1),1)*PARA.BGC.respFrac);
            
            
            %% Translocate mass to SOM pools :
            translocIn = zeros(nSoilCells,6);
            % From Lit 1, to SOM1
            translocIn(:, 4) = translocOut(:,1) + translocOut(:,2) + 0.93*translocOut(:,5) + translocOut(:,6);
            % From Lit 2 and SOM 1 to SOM 2:
            translocIn(:, 5) = translocOut(:, 3) + 0.9733*translocOut(:, 4);
            % From Lit3 and SOM2 to SOM 3:
            translocIn(:, 6) = 0.0267*translocOut(:, 4)+ 0.07*translocOut(:,5);
            
            % Calculate net mass change rate for each pool in each cell.
            Crate = Crate -decompRates +translocIn;
            
            
            % Calculate the rate of O2 consumption and CO2 production:
            %Calculate the respired fractions:
            respRate = decompRates*PARA.BGC.respFrac';
            
            O2rate = -1000*32/12 * sum(respRate,2);
            CO2rate = -44/32 * O2rate;
            
            gasRate(GRID.soil.cT_domain(GRID.BGC.cT_domain), 1:2) = gasRate(GRID.soil.cT_domain(GRID.BGC.cT_domain), 1:2) + [O2rate CO2rate];
            totalCO2rate = sum(CO2rate);
            
            if any(O2rate > 0)
                error('increasing O2')
            end
            GRID.BGC.kappa_O2 = kappaO2;
            GRID.BGC.kappa_CO2 = kappa_CO2;
            GRID.BGC.kappaWater = kappa_water;
        end
        
        function [ Crate, gasRate, CO2anaero, CH4anaero] = anaerobicRate( GRID, PARA, T, Crate, gasRate )
            %ANAEROBIC Calculates the rate of methane production and the
            %associated change in SOC.
            
            R = 8.314;   %[J / mol K] Universal gas constant
            %Find the subdomain of the BGC domain that consists of only soil cells
            decompDomain = GRID.soil.cT_domain(GRID.BGC.cT_domain);
            
            T = T(GRID.soil.cT_domain & GRID.BGC.cT_domain);   % Temperature in the soil cells in the BGC domain
            
            %Calculate the mass of O2 per porous volume (Grid cell volume cancels out):
            O2massPerPorVol = GRID.BGC.m_gases(decompDomain,1)./GRID.BGC.thetaSat(decompDomain);
            
            %Calculate the temperature scalar:
            Q10 = PARA.BGC.Q10_CH4gen.^((T-PARA.BGC.Tref_CH4gen)/10);
            % Q10(T<0) = 0;
            
            % cCO2cutOff = PARA.BGC.alpha_CO2 * 5.36*10^5 ./(T+273.15);          %[g/m^3]
            % cCO2 = GRID.BGC.c(decompDomain,2);
            %
            % kappa_CO2 = zeros(GRID.BGC.nSoilCells,1);
            % belowLimit = cCO2< cCO2cutOff;
            % if any(belowLimit)
            %     kappa_CO2(belowLimit) = (1-(cCO2(belowLimit)./cCO2cutOff(belowLimit)).^4);     %[-]
            % end
            kappa_CO2 = 1;
            
            %
            % Calculate the threshold pressure.
            P_threshold = GRID.BGC.waterP(decompDomain)*PARA.BGC.modelledGasFraction;
            
            %Calculate the partial pressures of all modelled gases:
            sumPartials = R*(T+273.15).*(GRID.BGC.c(decompDomain,:)*[32; 44; 16].^(-1));
            
            kappaPressure = max(1-(sumPartials./P_threshold).^4, 0);
            
            % Scale the rate down in proportion to the ice-filled porous volume:
            kappaIce = max(1- GRID.BGC.thetaIceSoil./GRID.BGC.thetaSatSoil,0);
            kappaO2 = exp(-O2massPerPorVol/PARA.BGC.O2star);
            
            
            % Calculate the rate scalar due to oxygen, temperature and pools:
            ratePerPool = (kappaO2.*Q10.*kappa_CO2.*kappaPressure.*kappaIce*PARA.BGC.f_pH_CH4gen)*(PARA.BGC.k_pool/5);
            decompRate  = GRID.BGC.C_mass.*ratePerPool;
            
            % Calculate the rate of change of mass of CH4 and CO2 [g m-2 s^-1
            CH4rate = 16/12 *PARA.BGC.molCH4perMolC *decompRate*transpose(PARA.BGC.respFrac)*1000;
            CO2rate = 44/16 *(1-PARA.BGC.molCH4perMolC)*decompRate*transpose(PARA.BGC.respFrac)*1000;
            
            
            %Calculate the masses that are translocated from each pool to another pool:
            translocOut =decompRate.*(1-ones(size(decompRate,1),1)*PARA.BGC.respFrac);
            
            
            % Translocate mass to SOM pools :
            translocIn = zeros(size(decompRate));
            % From Lit 1, to SOM1
            translocIn(:, 4) = translocOut(:,1) + translocOut(:,2) + 0.93*translocOut(:,5) + translocOut(:,6);
            % From Lit 2 and SOM 1 to SOM 2:
            translocIn(:, 5) = translocOut(:, 3) + 0.9733*translocOut(:, 4);
            % From Lit3 and SOM2 to SOM 3:
            translocIn(:, 6) = 0.0277*translocOut(:, 4)+ 0.07*translocOut(:,5);
            
            
            % Calculate the new masses:
            Crate = Crate - decompRate + translocIn;
            gasRate(decompDomain,2:3) = gasRate(decompDomain,2:3) + [CO2rate CH4rate ];
            CO2anaero = sum(CO2rate);
            CH4anaero = sum(CH4rate);
        end
        
        function [ gasRate, totCH4oxrate] = CH4Oxidation( GRID, PARA, T, gasRate)
            %CH4OXIDATION Calculates the rate of methane oxidation (methanotrophy).
            
            
            %Temperature in the the soil cells within the BGC domain:
            T = T(GRID.BGC.cT_domain & GRID.soil.cT_domain);
            
            % Vertical thickness of the soil cells.
            cellThicknesses = GRID.general.K_delta(GRID.BGC.cT_domain & GRID.soil.cT_domain);
            
            subdomain = GRID.soil.cT_domain(GRID.BGC.cT_domain);
            
            %Read the water-limit for the soil cells. This is calculated in the aerial decomposition module
            kappa_water = GRID.BGC.kappa_water;
            
            %Volume of air and water in each grid cell [m^3/m2 = m]
            totVolume = GRID.BGC.thetaSat(subdomain).*cellThicknesses;                %[m]
            %Calculate the overall concentrations in the pores:
            mPerPorV = GRID.BGC.m_gases(subdomain,:)./( totVolume*ones(1, GRID.BGC.ngases) );  %[g/m^3]
            
            %Read and calculate rate scalars:
            
            kappaO2 = mPerPorV(:,1)./(PARA.BGC.HT_halfRO2 + mPerPorV(:,1));       %[-]
            kappaCH4 = mPerPorV(:,3)./(PARA.BGC.halfR_CH4_CH4trophy+mPerPorV(:,3));
            
            Q10 = PARA.BGC.Q10_CH4trophy.^((T-PARA.BGC.Tref_CH4trophy)/10);
            
            % Potential mass of CH4 consumed during time step, if not O2-limited:
            CH4rate  = -PARA.BGC.Rmax_CH4trophy  .* kappaCH4 .* kappaO2 .* Q10 .* kappa_water;
            O2rate = 2*32/16*CH4rate;
            CO2rate = -O2rate*44/16;
            
            gasRate(subdomain,:) = gasRate(subdomain,:) +  [O2rate CO2rate CH4rate];
            totCH4oxrate = -sum(CH4rate);
            
            if any(O2rate > 0)
                error('increasing O2')
            end
        end
        
        
        function [ Crate ] = NPP(GRID, PARA, Crate )
            %NPPtoLitter allocates fresh litter to litter pools.
            
            totLit = PARA.BGC.annualNPP/31557600; %Litter flux rate, kg/s
            
            rootLit = totLit*PARA.BGC.subSurfNPPfrac;
            surfLit = max(0,totLit-rootLit);

            Crate(1,:) = Crate(1,:) + surfLit*PARA.BGC.NPP_f_litter;
            Crate(1:length(GRID.BGC.rootFrac),:) = Crate(1:length(GRID.BGC.rootFrac),:) + rootLit * (GRID.BGC.rootFrac*PARA.BGC.NPP_f_litter);
            
        end


        function rootFrac = rootDistribution( GRID, PARA )
            %ROOTDISTRIBUTION calculates the relative distribution of roots between
            %soil cells.
            dbstop if error
            
            ra = PARA.BGC.PFT.root_ra;
            rb = PARA.BGC.PFT.root_rb;
            
            soilz = GRID.general.K_grid(GRID.soil.K_domain);
            soilz = soilz - soilz(1);
            z = soilz(soilz <= PARA.BGC.rootDepth);
            
            rootFrac = zeros(size(z));
            if length(z) > 1
                rootFrac(1:end-1) = 0.5*( exp(-ra*z(1:end-1))+exp(-rb*z(1:end-1))-exp(-ra*z(2:end))-exp(-rb*z(2:end)) );
                rootFrac(end)= 0.5* (exp(-ra*z(end-1)+exp(-rb*z(end-1))));
            else
                rootFrac = 1;
            end
            % make it add up to 1:
            rootFrac = rootFrac+(1-sum(rootFrac))*rootFrac/sum(rootFrac);
        end

        function [ GRID ] = sedimentationSimpleNPP( GRID, PARA )  %must be handled in IA class
            
            % Makes grid cells expand of contract in response to SOM content and SOM density.
            
            GRID.BGC.bulkOrgThickness= GRID.BGC.C_mass/PARA.BGC.SOCperSOM * (PARA.BGC.poolBulkDens.^-1);
            
            %Update the inter-cellular boundaries below the soil surface:
            iGRID = GRID.BGC.K_domain_lb;
            iThickness = length(GRID.BGC.bulkOrgThickness);
            while iGRID > GRID.soil.K_domain_ub
                GRID.general.K_grid(iGRID-1) = GRID.general.K_grid(iGRID) - GRID.BGC.bulkOrgThickness(iThickness);
                iGRID = iGRID -1;
                iThickness = iThickness - 1;
            end
            
            % Update the inter-cellular boundaries above the soil surface:
            while iGRID > 1
                GRID.general.K_grid(iGRID-1) = GRID.general.K_grid(iGRID)-GRID.general.K_delta(iGRID-1);
                iGRID = iGRID-1;
            end
            
            % Recalculate the locations of the midpoints of the cells:
            GRID.general.cT_grid = 0.5*(GRID.general.K_grid(1:end-1) + GRID.general.K_grid(2:end));
            
            % Update the spacings:
            GRID.general.K_delta = GRID.general.K_grid(2:end)-GRID.general.K_grid(1:end-1);
            GRID.general.cT_delta = GRID.general.cT_grid(2:end) - GRID.general.cT_grid(1:end-1);
            
            % Recalculate the root distribution
            GRID.BGC.rootFrac = rootDistribution(GRID, PARA);
        end
        
        function [ newCrate , difMaxStep] = SOCdiffusionRate( GRID, PARA, Crate )
            % Read the depths of the cells from the lowest air cell to the lowest
            % BGC cell:
            z = GRID.general.cT_grid(GRID.soil.cT_domain_ub: GRID.BGC.cT_domain_lb);
            
            dz = z(2:end)-z(1:end-1);
            
            thicknesses = GRID.general.K_delta(GRID.soil.cT_domain & GRID.BGC.cT_domain);
            c = GRID.BGC.C_mass./(thicknesses*ones(1,size(GRID.BGC.C_mass,2)));   % kg/m3
            dc = c(2:end,:)-c(1:end-1,:); % kg/m3
            
            D = PARA.BGC.turbD; %Diffusivity   % m2 s-1
            
            FluxesDown = -D*dc./(dz*ones(1,size(GRID.BGC.C_mass,2)));  % m2/s * kg/m3/ m = kg/m2/s
            netInFlux = [-FluxesDown(1,:); FluxesDown(1:end-1,:)-FluxesDown(2:end,:); FluxesDown(end,:)]; % kg/m2/s
            
            newCrate = Crate +  netInFlux;  % kg/m2/s
            
            difMaxStep = 0.99*.5*min((dz.^2)/D);
        end

        
        
    end
end

