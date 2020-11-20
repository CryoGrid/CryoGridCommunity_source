%========================================================================
% CryoGrid TIER1 library class for functions related to soil gas

% S. Westermann, October 2020
%========================================================================

classdef SOIL_GAS < BASE
    
    methods
        
        function [ f ] = BubbleVolumeSolver(Vb, T, P_threshold, V_GC, c_row, epsilon_row )
            %not sure what this is for
            
            R = 8.314; %[J / mol K] Universal gas constant
            
            f = R*T*((c_row.*exp(-Vb*(epsilon_row*V_GC).^(-1)))*[32; 44; 16].^(-1))-P_threshold;
            %[J / mol K]*K(g/m3* exp(m^3*m^(-3)) * (g/mol)^-1) = Pa]
        end
        
        function [ GRID ] = cFromMass( GRID )
            %CFROMMASS Calculates the headspace concentration from the mass
            %concentration of each gas and the solubility-weighted porosity epsilon:
            %
            GRID.BGC.c = GRID.BGC.m_gases./(GRID.BGC.epsilon.*(GRID.general.K_delta(GRID.BGC.cT_domain)*ones(1,GRID.BGC.ngases)));
            % Avoid NaN when they can be excused:
            GRID.BGC.c(GRID.BGC.m_gases == 0 & GRID.BGC.epsilon == 0) = 0;
        end
        
        function [ GRID ] = massFrom_c( GRID )
            %CFROMMASS Calculates the headspace concentration from the mass
            %concentration of each gas and the solubility-weighted porosity epsilon:
            %
            GRID.BGC.m_gases = GRID.BGC.c.*(GRID.BGC.epsilon.*(GRID.general.K_delta(GRID.BGC.cT_domain)*ones(1,size(GRID.BGC.m_gases,2))));
        end

        
        
        function ground = soil_gas_diffusivities(ground)         %parameters for free-air diffusivities from LPj-WhyMe (Wania et al, 2010).
            %DIFFUSIVITIES This function calculates the diffusion
            %coefficients of the gases:
            T = ground.STATVAR.T;
            
            %Calculate the free-air and free-water diffusivities in the BGC domain,
            airD(:,1) = (.1875+0.0013*T)*10^-4;    %[m2/s] O2 diffusivity in free air
            airD(:,2) = (.1325+0.0009*T)*10^-4;    %[m2/s] CO2 diffusivity in free air
            airD(:,3) = (.1759+0.00117*T)*10^-4;   %[m2/s] CH4 diffusivity in free air
            
            ground.STATVAR.gas_diffusivity = airD.*(ground.STATVAR.air .* ones(1, size(airD,2)) + 10^-4 .* ground.STATVAR.water .* ones(1, size(airD,2)).* ground.STATVAR.Henry_sol);
            
%             % Set diffusivities in cells with significant amounts of ice to zero.
%             D((GRID.BGC.thetaAir+GRID.BGC.thetaWater)./GRID.BGC.thetaSat < 0.5,:) = 0;
%             %Calculate the effective diffusivity as minima of the diffusion coeffients
%             %in adjacent cells
%             % D = [min(D(1,:), D(2,:)); min( D(1:end-2,:), min( D(2:end-1,:), D(3:end,:))); min(D(end-1,:), D(end,:))];
% 
%             % Calculate the mean aerial diffusion coeffient for each gas in the cells
%             % containing aerenchyma:
%             nSnowCells = sum(double(GRID.snow.cT_domain));  %number of snow cells:
%             nArnchma  = min(sum(double(GRID.BGC.cT_domain)), length(GRID.BGC.rootFrac));   % number of soil cells with aerenchyma
%             ROI = nSnowCells+1:nSnowCells+nArnchma;
%             GRID.BGC.D_arnchma = max(0, sum(airD(ROI,3))/nArnchma);
%             
%             GRID.BGC.gasD = D;
%             if any( D < 0)
%                 error('Negative effective diffusivity')
%             end
            
        end
        
        
        function [ gasRate, difMaxStep ] = gasDifRateSurfEquil( GRID, gasRate, istest)
            % We use Fick's first law of diffusion for adjoining cells
            
            if ~istest
                D = GRID.BGC.gasD;
            else
                D0 = 1E-5;
                D = D0*ones(size(GRID.BGC.c));
            end
            
            z = GRID.general.cT_grid(GRID.BGC.cT_domain);
            
            
            % Set the diffusivity equal to the minimum diffusivity in adjoining pairs:
            D = min(D(1:end-1,:), D(2:end,:));
            
            % CAlculate concentration gradient
            dc = GRID.BGC.c(2:end,:) - GRID.BGC.c(1:end-1,:);
            dz = z(2:end)-z(1:end-1);
            dc_dz = dc./(dz*ones(1,3));
            
            % Upwards flux through each cell interface
            fluxInFromBelow = D.*dc_dz;   % m2/s * g/m3 / m = g/m2
            
            % The change in gas mass is the difference of flux in and flux out
            dm_dt = [ fluxInFromBelow(1,:); fluxInFromBelow(2:end,:)- fluxInFromBelow(1:end-1,:); -fluxInFromBelow(end,:)];
            
            % Update gasrate
            gasRate = gasRate +  dm_dt;
            
            % Set time step max
            difMaxStep =  .5*min(min((dz.^2)*ones(1,size(GRID.BGC.m_gases,2))./D));
        end
        
        
        function ground = Henry_solubility(ground, tile)
            %Henry_solubility This function calculates the dimmensionless Henry
            %solubility coefficients of O2, CO2 and CH4 given temperature T [K] and
            %pressure P [Pa].
            % Functions and coefficients derived from Sander, R. 2015. Compilation of
            % Henry's law constants (version 4.0) for water as solvent. Atmos. Chem.
            % Phys., 15, 4399-4981.

            P = ground.STATVAR.waterP; %atmospheric plus hydrostatic hydrostatic pressure, not sure
            
            P0 = 1.013E5;    %[kPa] Ambient pressure
            R  = ground.STATVAR.R; %8.314;      %[J K^-1 mol^-1] Universal gas constant
            T0 = 298.15;     %[K]   Standard temperature
            
            H0_O2  = 1.3E-5;    %[mol m^-3 Pa^-1]
            H0_CO2 = 3.5E-4;    %[mol m^-3 Pa^-1]
            H0_CH4 = 1.4E-5;    %[mol m^-3 Pa^-1]
            
            % T in the BGC domain
            T = ground.STATVAR.T;
            T = T + 273.15; %Convert from C to K
            
            Henry_sol(:,1) = H0_O2 *R*T.*P/P0.*exp(1500*(1./T-1/T0)); % O2
            Henry_sol(:,2) = H0_CO2*R*T.*P/P0.*exp(2400*(1./T-1/T0)); % CO2
            Henry_sol(:,3) = H0_CH4*R*T.*P/P0.*exp(1600*(1./T-1/T0)); % CH4
            
            ground.STATVAR.Henry_sol = Henry_sol;
            
        end
        
        function Cbefore  = totalCmass( GRID )
            %TOCMASS Calculates the total carbon mass in the GRID variable struct
            Cnow = 1000*sum(GRID.BGC.C_mass(:)) ...                  %SOC [g]
                + sum(GRID.BGC.m_gases(:,2:3)*[12/44; 12/16])...     %C in CO2 and CH4, respectively [g]
                + sum(GRID.BGC.gasInSoilIce(:,2:3)*[12/44; 12/16]);
            Clost = GRID.BGC.surfaceDiffusion(2:3)* [12/44; 12/16]  ...; %Cdiffused through the surface [g]
                + GRID.BGC.arnchmaFlux* 12/16 ...
                + GRID.BGC.surfaceEbullition(2:3)* [12/44; 12/16];
            Cgained = GRID.BGC.mGasAtmToSnow(2:3)*[12/44; 12/16] + GRID.BGC.NPP*1000;
            
            Cbefore = Cnow-Cgained+Clost;
        end
        
        function [ GRID ] = waterPressure( GRID, PARA, FORCING)  %treat this in the WATER classes
            
            %BGC_PRESSURE Calculates the sum of the atmospheric and hydrostatic
            %pressures. The water table height is considered as the height z_WT above
            %the parent material below which the soil is completely water-saturated.
            %   z:        Heigth above parent material                  (array)  [m]
            %   porosity: V-content of air, water, ice and nat.porosity [struct] [-]
            %   P_atm:    Atmospheric pressure                          (scalar) [Pa]
            
            %Declare constants and grids:
            water_density = 1000;    %[kg m^-3]
            g = 9.81;                %[m s^-2]
            
            P_atm = FORCING.i.p;    %atmospheric pressure
            
            %Depth constrained to the BGC domain:
            z = GRID.general.K_grid(GRID.BGC.K_domain);
            z = z(1:end-1);
            
            %Set the pressure equal to the atmospheric pressure and the variable for
            %water table depth empty:
            pressure = P_atm*ones(sum(double(GRID.BGC.cT_domain)),1);
            
            
            %Find the depth of the water table, z_WT:
            % isSaturated = abs(GRID.BGC.thetaSat(GRID.soil.cT_domain(GRID.BGC.cT_domain))-GRID.soil.cT_water(GRID.soil.cT_domain(GRID.BGC.cT_domain)))./...
            %     GRID.BGC.thetaSat(GRID.soil.cT_domain(GRID.BGC.cT_domain)) < PARA.BGC.satLimit;
            isSaturated = abs(GRID.BGC.thetaSat-GRID.BGC.thetaWater)./GRID.BGC.thetaSat < 1E-8;
            
            %If any of the grid cells is saturated:
            if any(isSaturated)
                waterTableIndex = find(isSaturated, 1); % Index of topmost saturated grid cell
                z_WT = z(waterTableIndex);         % Depth of water table
                %Add the hydrostatic pressure to the atmospheric pressure:
                pressure = pressure + water_density*g*(z-z_WT).*double(isSaturated);
            else
                waterTableIndex = NaN;
                z_WT = NaN;
            end
            % Generate output
            GRID.BGC.waterP = pressure;
            GRID.BGC.waterTableDepth = z_WT;
            GRID.BGC.waterTableIndex = waterTableIndex;
        end





    end
end

