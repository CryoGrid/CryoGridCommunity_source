function dD_dt = compact_windDrift_compiled(T, ice_snow, layerThick_snow, snow_mat1, wind_compaction_timescale, g)

%(ground.STATVAR.T(2:5,:), ground.STATVAR.ice_snow, ground.STATVAR.layerThick_snow, ground.TEMP.snow_mat1, tile.ENSEMBLE.STATVAR.wind_compaction_timescale, ground.CONST.g)



T = min(0,T);

eta_0 = 7.62237e6;
a_eta = 0.1;
b_eta = 0.023;
c_eta = 250;

rho_ice = 920;
rho_max = 350;

rho = ice_snow ./ max(1e-20, layerThick_snow) .* rho_ice;

%             rho(:,2)


stress = g .* rho_ice .* (cumsum(ice_snow - snow_mat1 .* ice_snow./2));

%             stress(:,2)


eta = eta_0 .*  rho ./ c_eta .* exp(-a_eta .* T + b_eta .* rho);

%             eta(:,2)

dD_dt =  - stress ./max(1e-10,eta) .* layerThick_snow; %compaction

%             dD_dt(:,2)


%dD_dt = dD_dt - PROFILE.snow_mat1 .* rho_ice .* PROFILE.D_ice_snow ./ max(1e-8, rho).^2.* (rho_max - min(rho , rho_max)) ./ tau_0 .* S_I;

dD_dt = dD_dt - snow_mat1 .* rho_ice .* ice_snow ./ max(1e-8, rho).^2.* (rho_max - min(rho, rho_max)) ./ repmat(wind_compaction_timescale .* 24.* 3600, 4,1);

%             dD_dt(:,2)


%dD_dt = dD_dt - repmat(double(sum(PROFILE.D_snow,1)>0.05),4,1) .* PROFILE.snow_mat1 .* rho_ice .* PROFILE.D_ice_snow ./ max(1e-8, rho).^2.* (rho_max - min(rho , rho_max)) ./PROFILE.wind_compaction_timescale;

