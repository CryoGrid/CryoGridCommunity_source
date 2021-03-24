%plots

plot(vegetation.mlcanopyinst.wind(1:vegetation.mlcanopyinst.ncan),height_ref(1:vegetation.mlcanopyinst.ncan))
plot(vegetation.mlcanopyinst.tleaf(2:vegetation.mlcanopyinst.ncan)-273.15,height_ref(2:vegetation.mlcanopyinst.ncan))
plot(vegetation.mlcanopyinst.tveg(1:vegetation.mlcanopyinst.ncan)-273.15,height_ref(1:vegetation.mlcanopyinst.ncan))
plot(vegetation.mlcanopyinst.tair(2:vegetation.mlcanopyinst.ncan)-273.15,height_ref(2:vegetation.mlcanopyinst.ncan))

plot(vegetation.mlcanopyinst.tleaf(2:vegetation.canopy.ntop)-273.15,height(2:vegetation.canopy.ntop))


%temperature profile vegetation
plot(vegetation.mlcanopyinst.tair(:,1:end)-273.15)
title({'Temperature of vegetation canopy layers';'Air = blue, Leafs = red, Canopy = yellow'})
xlabel('Canopy layer') 
ylabel('Temperature (°C)') 
hold
plot(vegetation.mlcanopyinst.tleaf(:,1:end,1)-273.15)
plot(vegetation.mlcanopyinst.tveg(:,1:end,1)-273.15)

%sunlit and shaded fractions of the canopy
plot(vegetation.flux.fracsun(:,:),height)
title('Fraction of sunny and shaded leaves')
xlabel('Fraction (0-1)') 
ylabel('Fraction of canopy height (1 = canopy top)')
hold
plot(vegetation.flux.fracsha(:,:),height)

%sum LAI (cumulative LAI)
plot(vegetation.canopy.sumlai(:,:))
title('Sum of LAI')
xlabel('Canopy layer') 
ylabel('LAI') 

%leaf temperature profile
plot(vegetation.mlcanopyinst.tleaf(:,:,1))
plot(vegetation.mlcanopyinst.tleaf(:,:,2))


plot(vegetation.mlcanopyinst.tair(:,:))

%Leaf absorbed solar radiation
plot(vegetation.flux.swleaf(:,1:end,2))
plot(vegetation.flux.swleaf(:,1:end,1))
title({'Leaf absorbed solar radiation';'VIS = blue, NIR = red'})
xlabel('Canopy layer') 
ylabel('Flux (W/m^2 leaf)') 

plot(vegetation.flux.lwleaf(:,:,1))

% vegetation.mlcanopyinst.sw_prof     = sw_prof      ; % Canopy layer absorbed solar radiation (W/m2)
% vegetation.mlcanopyinst.ir_prof     = ir_prof      ; % Canopy layer absorbed longwave radiation (W/m2)
% vegetation.mlcanopyinst.rn_prof     = rn_prof      ; % Canopy layer net radiation (W/m2)
% vegetation.mlcanopyinst.st_prof     = st_prof      ; % Canopy layer storage heat flux (W/m2)
% vegetation.mlcanopyinst.sh_prof     = sh_prof      ; % Canopy layer sensible heat flux (W/m2)
% vegetation.mlcanopyinst.lh_prof     = lh_prof      ; % Canopy layer latent heat flux (W/m2)
% vegetation.mlcanopyinst.et_prof     = et_prof      ; % Canopy layer water vapor flux (mol H2O/m2/s)
% vegetation.mlcanopyinst.fc_prof     = fc_prof      ; % Canopy layer CO2 flux (umol CO2/m2/s)

figure
plot(vegetation.mlcanopyinst.sh_prof(:,1:end))
title({'Canopy layer fluxes (W/m2)';'Sensible heat flux = blue, Storage heat flux = red'})
xlabel('Canopy layer') 
ylabel('Flux (W/m^2)') 
hold
plot(vegetation.mlcanopyinst.st_prof(:,1:end))

figure
plot(vegetation.mlcanopyinst.lh_prof(:,1:2))
title({'Canopy layer fluxes (W/m2)'; 'Latent heat flux = blue, Net radiation = red'})
xlabel('Canopy layer') 
ylabel('Flux (W/m^2)') 
hold
plot(vegetation.mlcanopyinst.rn_prof(:,1:end))

plot(vegetation.mlcanopyinst.tleaf(:,1:end,1)-273.15)
plot(vegetation.mlcanopyinst.tleaf_old(:,1:end,1)-273.15)

plot(vegetation.mlcanopyinst.rnleaf(:,1:end,1))
