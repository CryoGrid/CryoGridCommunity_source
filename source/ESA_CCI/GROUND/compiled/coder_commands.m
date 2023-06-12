codegen advance_E_compiled -args {zeros(44,2100), zeros(44,2100),0}
codegen compact_windDrift_compiled -args {zeros(4,2100), zeros(4,2100), zeros(4,2100), zeros(4,2100), zeros(1,2100), 0}
codegen dE_dt_compiled -args {zeros(44,2100), zeros(44,2100), zeros(45,2100), zeros(44,2100)}
codegen get_Tground_compiled -args {zeros(44,2100), zeros(41,2100), zeros(41,2100), zeros(41,2100), zeros(41,2100)}
codegen thermCond_compiled -args {zeros(45,2100), zeros(41,2100), zeros(41,2100), zeros(41,2100), zeros(41,2100)}
codegen thermCond_eff_compiled -args {zeros(45,2100), zeros(45,2100)}
codegen thermCond_snow_compiled -args {zeros(4,2100), zeros(4,2100)}

%not successful
% codegen advance_E_compiled -args {single(zeros(44,2100)), single(zeros(44,2100)),single(0)}
% codegen dE_dt_compiled -args {single(zeros(44,2100)), single(zeros(44,2100)), single(zeros(45,2100)), single(zeros(44,2100))}
% codegen get_Tground_compiled -args {single(zeros(44,2100)), single(zeros(41,2100)), single(zeros(41,2100)), single(zeros(41,2100)), single(zeros(41,2100))}
% codegen thermCond_compiled -args {single(zeros(45,2100)), single(zeros(41,2100)), single(zeros(41,2100)), single(zeros(41,2100)), single(zeros(41,2100))}
% codegen thermCond_eff_compiled -args {single(zeros(45,2100)), single(zeros(45,2100))}

codegen get_boundary_condition_u_compiled -args {ground.STATVAR.T, tile.FORCING.TEMP.surfT, tile.ENSEMBLE.PARA.ensemble_size, tile.FORCING.TEMP.melt_bare, tile.ENSEMBLE.STATVAR.melt_fraction, tile.FORCING.TEMP.melt_forest, ground.CONST.day_sec, tile.timestep, tile.FORCING.TEMP.snowfall, tile.ENSEMBLE.STATVAR.snowfall_factor, ground.STATVAR.ice_snow, ground.STATVAR.upper_cell, ground.TEMP.d_energy,ground.TEMP.snow_mat1, ground.CONST.L_f, ground.CONST.c_i, ground.STATVAR.layerThick_snow, tile.ENSEMBLE.STATVAR.wind_speed_class}
codegen get_boundary_condition_u_compiled -args {zeros(45,2100), zeros(1,2100), zeros(1,2100), 0,0, zeros(1,2100), zeros(4,2100), zeros(1,2100), zeros(44,2100), zeros(4,2100), 0, 0, zeros(4,2100), zeros(1,2100)}


get_boundary_condition_u_compiled(T, surfT, melt, ...
 day_sec, timestep, snow_fall, ice_snow, upper_cell, d_energy, snow_mat1, L_f, c_i, layerThick_snow, wind_speed_class)