function [ground] = initialize_pftcon(ground)

ground.STATVAR.vegetation.pftcon.beta_neutral_max = 0.35;
ground.STATVAR.vegetation.pftcon.cr = 0.3;
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% 1. CanopyNitrogen
% Input
ground.STATVAR.vegetation.pftcon.vcmaxpft = 43;       %= pftcon.vcmaxpft;           % Maximum carboxylation rate at 25C (umol/m2/s)
ground.STATVAR.vegetation.pftcon.c3psn = 1;           %= pftcon.c3psn    ;          % Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant
ground.STATVAR.vegetation.pftcon.emleaf = 0.98;       % Leaf emissivity

%http://www.cgd.ucar.edu/tss/clm/pfts/pft-physiology.htm
%https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/pftdata/pft-physiology-cn.c070212.readme
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% 1. CanopyWater
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%CiFunc (not translated yet)
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%CiFuncGS (not translated yet)
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Forest Radiation
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%FrictionVelocity
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%ft
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%fth
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%fth25
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%GetPsiCLM
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%LatVap
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Leaf Boundry Layer DONE
%Input
ground.STATVAR.vegetation.pftcon.dleaf = 0.04;  %http://www.cgd.ucar.edu/tss/clm/pfts/pft-physiology.htm   % Leaf dimension (m)

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Leaf Fluxes DONE
%1.1 Leaf fluxes
%1.2 T leaf function
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Leaf Temperature DONE
%1.1 Leaf Temperature and Energy Fluxes
%1.2 Leaf Heat Capacity
%Input
ground.STATVAR.vegetation.pftcon.slatop = 0.008;              % Bonan 2018, Table 1    % Specific leaf area at top of canopy (m2/gC)
%-----------------------------------------------------------------------
%LeafWaterPotential
% ground.STATVAR.vegetation.pftcon.capac = ones(1);     % Plant capacitance (mmol H2O/m2 leaf area/MPa)

%LeafHeatCapacity
ground.STATVAR.vegetation.pftcon.cpliq = 4188.;                        %SHR_CONST_CPFW specific heat of fresh h2o ~ J/kg/K %https://github.com/ESCOMP/ctsm/blob/84a1ed38ea4d8c89d6f51bcd8b377fb49b8311f6/tools/mksurfdata_map/src/shr_const_mod.F90
ground.STATVAR.vegetation.pftcon.cpbio = 1396.;                           % 800-1480 J/kg/C (specific heat of dry-wet soil, https://www.engineeringtoolbox.com/specific-heat-capacity-d_391.html) %use clm_varcon
ground.STATVAR.vegetation.pftcon.fcarbon = 0.5;                          % Fraction of dry biomass that is carbon
ground.STATVAR.vegetation.pftcon.fwater = 0.7;                           % Fraction of fresh biomass that is water
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Monin ObukIniCLM DONE
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Obukhov DONE
%Input
ground.STATVAR.vegetation.pftcon.z0mr = 0.055;        %https://ldas.gsfc.nasa.gov/gldas/data/GLDAS1_vegparam_tbl_clm2l.pdf  % Ratio of momentum roughness length to canopy top height
ground.STATVAR.vegetation.pftcon.displar = 0.67;      %https://ldas.gsfc.nasa.gov/gldas/data/GLDAS1_vegparam_tbl_clm2l.pdf  % Ratio of displacement height to canopy top height

%--------------------------------------------------------------------------------------------------------------------- --------------------------------------------------------------
%Photosynthesis DONE
%1.1 Leaf-level photosynthesis variables

%1.1 PhotosynthesisParam (p, mlcanopy_inst)

%1.2 Leaf photosynthesis and stomatal conductance
%Input
%c3psn = 1 ;                    % Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant
ground.STATVAR.vegetation.pftcon.g1opt  = 6;  % Slope of conductance (mp) https://ldas.gsfc.nasa.gov/gldas/data/GLDAS1_vegparam_tbl_clm2l.pdf  % Ball-Berry slope of conductance-photosynthesis relationship, unstressed
ground.STATVAR.vegetation.pftcon.g0opt  = 0.01;  % Ball-Berry minimum leaf conductance, unstressed (mol H2O/m2/s)    Bonan Matlab leaf physiology params 12_01

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%PlantResistance
%   gplant    = pftcon.gplant       ;   % Stem (xylem-to-leaf) hydraulic conductance (mmol H2O/m2 leaf area/s/MPa)
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%SoilResistance
% Bonan book p. 221 (resistivity) and p. 35
ground.STATVAR.vegetation.pftcon.root_radius = 0.29e-03;          % Fine root radius (m)
ground.STATVAR.vegetation.pftcon.root_density = 0.31e06;           % Fine root density (g biomass / m3 root)     % http://www.cesm.ucar.edu/events/wg-meetings/2015/presentations/bgcwg+lmwg/drewniak.pdf
ground.STATVAR.vegetation.pftcon.root_resist = 25; %200000;         % Hydraulic resistivity of root tissue (MPa.s.g/mmol H2O) rootresist(MPa s g mol-1(biomass based)) 30-150 Optimal value % https://yncenter.sites.yale.edu/sites/default/files/wei_nan_may_2018.pdf
end