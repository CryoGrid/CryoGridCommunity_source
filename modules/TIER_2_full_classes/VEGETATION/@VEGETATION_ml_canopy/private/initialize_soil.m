% --- Soil variables
function [ground] = initialize_soil(ground)

ground.STATVAR.vegetation.soilvar.nsoi = 10; %same as layers in the GROUND class
ground.STATVAR.vegetation.soilvar.rootfr  = [0; 0.25./3.*ones(3,1); 0.75./6 .* ones(6,1)]'; %Simone, change root fraction here!
ground.STATVAR.vegetation.soilvar.dz = zeros(1, ground.STATVAR.vegetation.soilvar.nsoi);
ground.STATVAR.vegetation.soilvar.h2osoi_vol = zeros(1, ground.STATVAR.vegetation.soilvar.nsoi);
ground.STATVAR.vegetation.soilvar.h2osoi_ice = zeros(1, ground.STATVAR.vegetation.soilvar.nsoi);

% values from here: https://slideplayer.com/slide/12736113/76/images/21/Clapp+and+Hornberger+(1978)+parameters+for+Soil+Moisture+Characteristic+functions+based+on+analysis+of+1845+soils.+Values+in+parentheses+are+standard+deviations..jpg
% ground.STATVAR.vegetation.soilvar.dz           = 0.4;  %0.0451; %0.0175;   %ground.STATVAR.vegetation.soilvar.dz         ;  % Soil layer thickness (m)
% ground.STATVAR.vegetation.soilvar.watsat       = 0.2;  %ground.STATVAR.vegetation.soilvar.watsat     ;  % Soil layer volumetric water content at saturation (porosity)

% ground.STATVAR.vegetation.soilvar.transp_per_layer = [0.6,0.3,0.1];
% ground.STATVAR.vegetation.soilvar.hksat        = [1.056,0.938,0.0147]; %3; %0.00001; %0.01; %3.;   %soilstate_inst.hksat_col      ;  % Soil layer hydraulic conductivity at saturation (mm H2O/s)
% ground.STATVAR.vegetation.soilvar.bsw          = [4.05,4.38,10.40];   %soilstate_inst.bsw_col        ;  % Soil layer Clapp and Hornberger "b" parameter

% % !!!!!!!!!!!!!!!!!!!!!!!!! SIMONE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ground.STATVAR.vegetation.soilvar.rootfr       = [0.5,0.25,0.25];  % !!!!!!!!!!!!!!!!!!!!!!!!! SIMONE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %soilstate_inst.rootfr_patch   ;  % Fraction of roots in each soil layer
% ground.STATVAR.vegetation.soilvar.h2osoi_vol   = [0.2,0.2,0.2];  %waterstate_inst.h2osoi_vol_col;  % Soil layer volumetric water content (m3/m3)
% ground.STATVAR.vegetation.soilvar.h2osoi_ice   = [0.,0.,0.];   %waterstate_inst.h2osoi_ice_col;  % Soil layer ice lens (kg/m2)

% ground.STATVAR.vegetation.soilvar.t_top_surfacecell = ground.STATVAR.vegetation.mlcanopyinst.tref+10;
% ground.STATVAR.vegetation.soilvar.dz_topsurfacecell = 0.05;
% ground.STATVAR.vegetation.soilvar.thk_topsurfacecell = 2.0;
% 
% % % Soil layer thickness (m) -> Thickness of each layer
% % ground.STATVAR.vegetation.soilvar.dz = [0.15, 0.3, 0.5];
% 
% ground.STATVAR.vegetation.soilvar.z = [0.1, 0.2, 0.7];     % Soil layer depth (m)
% ground.STATVAR.vegetation.soilvar.zi = [0.05, 0.15, 0.45];     % Soil layer depth at layer interface (m)
% ground.STATVAR.vegetation.soilvar.nsoi = 10;
% 
% % Soil layer thickness (m)
% ground.STATVAR.vegetation.soilvar.dz = [0.1, 0.1, 0.5];
% % 
% % Simone: why was this in snow module?
% ground.STATVAR.vegetation.soilvar.t_soisno  = [ground.STATVAR.vegetation.mlcanopyinst.tref+10,ground.STATVAR.vegetation.mlcanopyinst.tref+10,ground.STATVAR.vegetation.mlcanopyinst.tref+10]; %273.15+20; %273.15+15;  % Soil temperature (K)                    %260 woher auch immer??? 
% ground.STATVAR.vegetation.soilvar.thk = [2.0,2.0,2.0]; %-> Now exchanged in IA class
% 
% % Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface)
% ground.STATVAR.vegetation.soilvar.z_plus_onehalf(1) = -ground.STATVAR.vegetation.soilvar.dz(1);
% for i = 2:ground.STATVAR.vegetation.soilvar.nsoi
%     ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i) = ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i-1) - ground.STATVAR.vegetation.soilvar.dz(i);
% end
% 
% % Soil depth (m) at center of layer i (negative distance from surface)
% 
% ground.STATVAR.vegetation.soilvar.z(1) = 0.5 * ground.STATVAR.vegetation.soilvar.z_plus_onehalf(1);
% for i = 2:ground.STATVAR.vegetation.soilvar.nsoi
%     ground.STATVAR.vegetation.soilvar.z(i) = 0.5 * (ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i-1) + ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i));
% end
% %
% % %Thickness between between z(i) and z(i+1)
% 
% for i = 1:ground.STATVAR.vegetation.soilvar.nsoi-1
%     ground.STATVAR.vegetation.soilvar.dz_plus_onehalf(i) = ground.STATVAR.vegetation.soilvar.z(i) - ground.STATVAR.vegetation.soilvar.z(i+1);
% end
% ground.STATVAR.vegetation.soilvar.dz_plus_onehalf(ground.STATVAR.vegetation.soilvar.nsoi) = 0.5 * ground.STATVAR.vegetation.soilvar.dz(ground.STATVAR.vegetation.soilvar.nsoi);
% % --- Initial conditions
% 
% % Initial soil temperature (K) and unfrozen and frozen water (kg H2O/m2)
% 
% for i = 1:ground.STATVAR.vegetation.soilvar.nsoi
%     
%     % Temperature
%     tmean = ground.STATVAR.vegetation.mlcanopyinst.tref-ground.STATVAR.vegetation.physcon.tfrz; %15.0;             % Mean daily air temperature (C)
%     ground.STATVAR.vegetation.soilvar.tsoi(i) = tmean + ground.STATVAR.vegetation.physcon.tfrz;
%     
%     % Soil water at saturation (kg H2O/m2)
%     
%     % Watsat and soil_texture taken from sp_05_01.m
%     % Volumetric soil water content at saturation (porosity)
%     % (Clapp and Hornberger. 1978. Water Resources Research 14:601-604)
%     
%     ground.STATVAR.vegetation.soilvar.watsat = [0.8, 0.8, 0.6]; % simone: Watsat ist die Porosität, die muss zumindest am Anfang mit der von CryoGrid übereinstimmen. Am besten würde die Porosität aus CryoGrid in die Vegetation eingelesen werden sollen. 
%     
%     h2osoi_sat = ground.STATVAR.vegetation.soilvar.watsat(i) * 1000 * ground.STATVAR.vegetation.soilvar.dz(i);
%     
%     % ground.STATVAR.vegetation.soilvar.soil_water_matric_potential = -(1/alpha)*(((0w - 0wr) / (0ws - 0wr))^(n /(1-n))-1)^(1/n);
%     
%     %soil parameters sandy clay loam (p.120, bonan book)
%     ground.STATVAR.vegetation.soilvar.alpha = 5.9; %cm -> m  0.059;
%     ground.STATVAR.vegetation.soilvar.n = 1.48;
%     ground.STATVAR.vegetation.soilvar.Owr(i) = ground.STATVAR.vegetation.soilvar.watsat(i);
%     ground.STATVAR.vegetation.soilvar.Ow(i) = ground.STATVAR.vegetation.soilvar.h2osoi_vol(i);
%     ground.STATVAR.vegetation.soilvar.Ows = 0;
%     
%     ground.STATVAR.vegetation.soilvar.soil_water_matric_potential(i) = -(1/ground.STATVAR.vegetation.soilvar.alpha)*(((ground.STATVAR.vegetation.soilvar.Ow - ground.STATVAR.vegetation.soilvar.Owr) / (ground.STATVAR.vegetation.soilvar.Ows - ground.STATVAR.vegetation.soilvar.Owr))^(ground.STATVAR.vegetation.soilvar.n /(1-ground.STATVAR.vegetation.soilvar.n)-1))^(1/ground.STATVAR.vegetation.soilvar.n); %m
%     
%     % Actual water content is some fraction of saturation. These are only used for soil
%     % thermal properties and phase change. Note the inconsistency with the use of soil
%     % water in the bucket model to calculate the soil wetness factor.
%     
%     satfrac = 0.4;
%     if (ground.STATVAR.vegetation.soilvar.tsoi(i) >= 273.15)
%         ground.STATVAR.vegetation.soilvar.h2osoi_ice(i) = 0;
%         ground.STATVAR.vegetation.soilvar.h2osoi_liq(i) = satfrac * h2osoi_sat;
%     else
%         ground.STATVAR.vegetation.soilvar.h2osoi_liq(i) = 0;
%         ground.STATVAR.vegetation.soilvar.h2osoi_ice(i) = satfrac * h2osoi_sat;
%     end
    
end


