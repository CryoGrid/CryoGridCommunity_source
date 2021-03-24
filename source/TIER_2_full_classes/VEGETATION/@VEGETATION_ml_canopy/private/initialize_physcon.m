
function [ground] = initialize_physcon(ground)

% https://escomp.github.io/ctsm-docs/doc/build/html/tech_note/Ecosystem/CLM50_Tech_Note_Ecosystem.html#table-physical-constants
%clm_arguments structure/container
% from sp_16_01
ground.STATVAR.vegetation.physcon.vkc = 0.4;                   % von Karman constant
ground.STATVAR.vegetation.physcon.grav = 9.80616;              % Gravitational acceleration (m/s2)
ground.STATVAR.vegetation.physcon.tfrz = 273.15;               % Freezing point of water (K)
ground.STATVAR.vegetation.physcon.hvap = 2.501e6;              % Latent heat of evaporation (J/kg)   = 45.06kJ/mol-1
ground.STATVAR.vegetation.physcon.mmh2o = 18.016 ./ 1000;       % Molecular mass of water (kg/mol)
ground.STATVAR.vegetation.physcon.mmdry = 28.966 ./ 1000;       % Molecular mass of dry air (kg/mol)
ground.STATVAR.vegetation.physcon.cpd = 1004.64;             % Specific heat of dry air at constant pressure (J/kg/K)
ground.STATVAR.vegetation.physcon.denh20  = 1000.;           % Density of fresh water ~ kg/m^3 %https://github.com/ESCOMP/ctsm/blob/84a1ed38ea4d8c89d6f51bcd8b377fb49b8311f6/tools/mksurfdata_map/src/shr_const_mod.F90 (https://github.com/ESCOMP/ctsm/blob/f2839e2f7c3374ee357508021594d13d5aa5534b/src/main/clm_varcon.F90)
ground.STATVAR.vegetation.physcon.rhoice = 917.;                % Density of ice (kg/m3)
ground.STATVAR.vegetation.physcon.hfus = 334000; %3.337e5;              % Heat of fusion for water at 0 C (J/kg)

% from sp_07_01
ground.STATVAR.vegetation.physcon.sigma = 5.67e-08;            % Stefan-Boltzmann constant (W/m2/K4)
% ground.STATVAR.vegetation.physcon.cwat = 4188.;              % Specific heat of water (J/kg/K)
% ground.STATVAR.vegetation.physcon.cice = 2117.27;            % Specific heat ice (J/kg/K)
ground.STATVAR.vegetation.physcon.tkwat = 0.57;                % Thermal conductivity of water (W/m/K)
ground.STATVAR.vegetation.physcon.tkice = 2.29;                % Thermal conductivity of ice (W/m/K)
ground.STATVAR.vegetation.physcon.hsub = ground.STATVAR.vegetation.physcon.hfus + ground.STATVAR.vegetation.physcon.hvap;     % Latent heat of sublimation (J/kg)

ground.STATVAR.vegetation.physcon.pi = 3.14159;

% from sp_12_02
ground.STATVAR.vegetation.physcon.visc0 = 13.3e-06;             % Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
ground.STATVAR.vegetation.physcon.Dh0 = 18.9e-06;               % Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
ground.STATVAR.vegetation.physcon.Dv0 = 21.8e-06;               % Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
ground.STATVAR.vegetation.physcon.Dc0 = 13.8e-06;               % Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)
ground.STATVAR.vegetation.physcon.rhowat = 1000;                % Density of water (kg/m3)
ground.STATVAR.vegetation.physcon.cpw = 1810; %46.;                   % Specific heat of water vapor at constant pressure (J/kg/K)
% ground.STATVAR.vegetation.physcon.cvwat = ground.STATVAR.vegetation.physcon.cwat * ground.STATVAR.vegetation.physcon.rhowat;  % Heat capacity of water (J/m3/K)
% ground.STATVAR.vegetation.physcon.cvice = ground.STATVAR.vegetation.physcon.cice * ground.STATVAR.vegetation.physcon.rhoice;  % Heat capacity of ice (J/m3/K)
ground.STATVAR.vegetation.physcon.rgasc =  8.3145; %8.314462618;    %[J/K/mole]           %8.314462618e3;        %!universal gas constant [J/K/kmole]

%Arguments and Uses

% VARCON variables from
%       https://github.com/ESCOMP/ctsm/blob/f2839e2f7c3374ee357508021594d13d5aa5534b/src/main/clm_varcon.F90
%       (clm_varcon.F90)
%       and
%       https://github.com/ESCOMP/ctsm/blob/84a1ed38ea4d8c89d6f51bcd8b377fb49b8311f6/tools/mksurfdata_map/src/shr_const_mod.F90
%       (shr_const_mod.F90)

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% 1. CanopyNitrogen DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% 1. CanopyWater DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------                    
%CiFunc (not translated yet)
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------                    
%CiFuncGS (not translated yet)
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%Forest Radiation DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%FrictionVelocity DONE
% PFT table Bonan
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%ft DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%fth DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%fth25 DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% %GetPsiCLM DONE
%             ground.STATVAR.vegetation.physcon.zetam = -1.574;         %CLM: transition point of flux-gradient relation (wind profile)
%             ground.STATVAR.vegetation.physcon.zetah = -0.465;         %CLM: transition point of flux-gradient relation (temperature profile)        
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%LatVap
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Leaf Boundry Layer DONE
%             ground.STATVAR.vegetation.physcon.b1 = 1.5;
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Leaf Fluxes DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Leaf Temperature DONE
%-----------------------------------------------------------------------
%LeafWaterPotential            
%-----------------------------------------------------------------------
%LeafWaterPotential
%             ground.STATVAR.vegetation.physcon.head = ground.STATVAR.vegetation.physcon.denh20*ground.STATVAR.vegetation.physcon.grav*1.e-06;  % Head of pressure  (MPa/m)
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Monin ObukIniCLM DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Obukhov DONE       
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Photosynthesis DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
%PlantResistance
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------            
%SoilResistance
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Satvap DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Scalar Profile
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%StabilityFunc1 DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%StabilityFunc2 DONE 
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%StomataEfficiency         
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%StomataFluxes         
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%StomataOptimization            
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
end