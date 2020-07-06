# CryoVeg
Created 20191001, CryoGrid branch
This is the modular version of CryoGrid which is based on an object-oriented programming paradigm, coupled to a multilayer canopy model (Bonan 2018)

20191001 Simone Stünzi

What has changed? 
- modules/INTERACTION/IA_HEAT_VEGETATION -> New interaction class
- modules/INTERACTION/get_IA_class -> New definition for interactions with vegetation module
- modules/@GROUND_vegetation -> New vegetation module based on https://github.com/gbonan/CLM-ml_v0
- modules/FORCING/@FORCING_seb/Nyurba_CCSM_rcp8_5_1979_2100.mat (not both necessary) -> Forcing file for the Nyurba site
- results/test_excel/test_excel.xlsx -> New class (GROUND_vegetation) and new stratigraphy

Model set-up and parameters used described in https://www.biogeosciences-discuss.net/bg-2020-201/#discussion

