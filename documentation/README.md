# CryoGrid Vegetation
Created 20191001, CryoGrid branch
This is the modular version of CryoGrid which is based on an object-oriented programming paradigm, coupled to a multilayer canopy model developed by Bonan et al., 2018

Bonan, G. B., Patton, E. G., Harman, I. N., Oleson, K. W., Finnigan, J. J., Lu, Y., and Burakowski, E. A.: Modeling canopy-induced turbulencein the Earth system: A unified parameterization of turbulent exchange within plant canopies and the roughness sublayer (CLM-ml v0),Geoscientific Model Development, 11, 1467–1496, https://doi.org/10.5194/gmd-11-1467-2018, 2018

20191001 Simone Stünzi

What has changed? 
- modules/INTERACTION/IA_HEAT_VEGETATION -> New interaction class
- modules/INTERACTION/get_IA_class -> New definition for interactions with vegetation module
- modules/@GROUND_vegetation -> New vegetation module based on https://github.com/gbonan/CLM-ml_v0
- modules/@GROUND_vegetation/private -> Translated functions from https://github.com/gbonan/CLM-ml_v0
- modules/@GROUND_vegetation/private/canopy_fluxes_multilayer -> Computed fluxes for each layer (sunlit and shaded)
- modules/FORCING/@FORCING_seb/Nyurba_CCSM_rcp8_5_1979_2100.mat -> Forcing file for the Nyurba site (N 63.18946°, E 118.19596°)
- results/test_excel/test_excel.xlsx -> New class (GROUND_vegetation) and new stratigraphy

Model set-up and parameters used described in https://www.biogeosciences-discuss.net/bg-2020-201/#discussion

