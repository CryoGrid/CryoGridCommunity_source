# CryoVeg
Created 20191001, private CryoGrid branch
This is the modular version of CryoGrid which is based on an object-oriented programming paradigm, coupled to a multilayer canopy model

20191001 Simone Stünzi

What has changed? 
- modules/INTERACTION/IA_HEAT_VEGETATION -> New interaction class
- modules/INTERACTION/get_IA_class -> New definition for interactions with vegetation module
- modules/@GROUND_vegetation -> New vegetation module
- modules/FORCING/@FORCING_seb/CG3_CCLM_forcing_90_101.mat / forcing/CG3_CCLM_forcing_90_101.mat (not both necessary) -> New forcing
- results/test_excel/test_excel.xlsx -> New class (GROUND_vegetation) and new stratigraphy


% For now I changed the start time (l.48) in the main file! 
