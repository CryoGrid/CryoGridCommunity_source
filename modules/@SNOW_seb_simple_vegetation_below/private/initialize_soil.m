% --- Soil variables
function [ground] = initialize_soil(ground)

% values from here: https://slideplayer.com/slide/12736113/76/images/21/Clapp+and+Hornberger+(1978)+parameters+for+Soil+Moisture+Characteristic+functions+based+on+analysis+of+1845+soils.+Values+in+parentheses+are+standard+deviations..jpg
% ground.STATVAR.vegetation.soilvar.dz           = 0.4;  %0.0451; %0.0175;   %ground.STATVAR.vegetation.soilvar.dz         ;  % Soil layer thickness (m)
% ground.STATVAR.vegetation.soilvar.watsat       = 0.2;  %ground.STATVAR.vegetation.soilvar.watsat     ;  % Soil layer volumetric water content at saturation (porosity)


ground.STATVAR.vegetation.soilvar.hksat        = 0.176; %3; %0.00001; %0.01; %3.;   %soilstate_inst.hksat_col      ;  % Soil layer hydraulic conductivity at saturation (mm H2O/s)
ground.STATVAR.vegetation.soilvar.bsw          = 4.05;   %soilstate_inst.bsw_col        ;  % Soil layer Clapp and Hornberger "b" parameter


% ground.STATVAR.vegetation.soilvar.smp_l        = 1.;   %soilstate_inst.smp_l_col      ;  % Soil layer matric potential (mm)
ground.STATVAR.vegetation.soilvar.rootfr       = 0.2;  %soilstate_inst.rootfr_patch   ;  % Fraction of roots in each soil layer

%volumetric water content sandy clay loam (p.120, bonan book)
ground.STATVAR.vegetation.soilvar.h2osoi_vol   = 0.10;  %waterstate_inst.h2osoi_vol_col;  % Soil layer volumetric water content (m3/m3)
ground.STATVAR.vegetation.soilvar.h2osoi_ice   = 0.;   %waterstate_inst.h2osoi_ice_col;  % Soil layer ice lens (kg/m2)
            
ground.STATVAR.vegetation.soilvar.thk = 1.2;%0.08;          % W/m*K     1.8;    % Thermal conductivity of the first snow/soil layer; Larch forest = https://lab.agr.hokudai.ac.jp/env/ctc_siberia/bibliography/pdf/Brouchkov2005.pdf
% ground.STATVAR.vegetation.soilvar.cv = 1180.;%1180.;           % soil heat capacity (J/m2/K) Soil heat capacity from CryoGrid3
            
ground.STATVAR.vegetation.soilvar.t_soisno  = ground.STATVAR.vegetation.mlcanopyinst.tref; %273.15+20; %273.15+15;  % Soil temperature (K)                    %260 woher auch immer??? 
             
% ground.STATVAR.vegetation.soilvar.smp_l = 0.15;        % Soil layer matric potential (mm)
% ground.STATVAR.vegetation.soilvar.snl = 1;           % Number of snow layers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ground.STATVAR.vegetation.soilvar.z = 0.05; %0.6; %0.20;            % Soil layer depth (m)
ground.STATVAR.vegetation.soilvar.zi = 0.0;%.14;            % Soil layer depth at layer interface (m)


% Number of layers in soil profile
% set in SetUpCanopy
% ground.STATVAR.vegetation.soilvar.nsoi = 1;

% Soil layer thickness (m)7.1006354171935350E-0080.03
ground.STATVAR.vegetation.soilvar.dz = [0.05, 0.05, 0.05, 0.05, 0.05]; %005, 0.25, 0.5]; %0175; %, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058];

% Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface)

ground.STATVAR.vegetation.soilvar.z_plus_onehalf(1) = -ground.STATVAR.vegetation.soilvar.dz(1);
for i = 2:ground.STATVAR.vegetation.soilvar.nsoi
    ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i) = ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i-1) - ground.STATVAR.vegetation.soilvar.dz(i);
end

% Soil depth (m) at center of layer i (negative distance from surface)

% ground.STATVAR.vegetation.soilvar.z(1) = 0.5 * ground.STATVAR.vegetation.soilvar.z_plus_onehalf(1);
% for i = 2:ground.STATVAR.vegetation.soilvar.nsoi
%     ground.STATVAR.vegetation.soilvar.z(i) = 0.5 * (ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i-1) + ground.STATVAR.vegetation.soilvar.z_plus_onehalf(i));
% end
% 
% %Thickness between between z(i) and z(i+1)
% 
% for i = 1:ground.STATVAR.vegetation.soilvar.nsoi-1
%     ground.STATVAR.vegetation.soilvar.dz_plus_onehalf(i) = ground.STATVAR.vegetation.soilvar.z(i) - ground.STATVAR.vegetation.soilvar.z(i+1);
% end
% ground.STATVAR.vegetation.soilvar.dz_plus_onehalf(ground.STATVAR.vegetation.soilvar.nsoi) = 0.5 * ground.STATVAR.vegetation.soilvar.dz(ground.STATVAR.vegetation.soilvar.nsoi);

% --- Initial conditions

% Initial soil temperature (K) and unfrozen and frozen water (kg H2O/m2)

for i = 1:ground.STATVAR.vegetation.soilvar.nsoi
    
    % Temperature
    tmean = ground.STATVAR.vegetation.mlcanopyinst.tref-ground.STATVAR.vegetation.physcon.tfrz; %15.0;             % Mean daily air temperature (C)
    ground.STATVAR.vegetation.soilvar.tsoi(i) = tmean + ground.STATVAR.vegetation.physcon.tfrz;
    
    % Soil water at saturation (kg H2O/m2)
    
    % Watsat and soil_texture taken from sp_05_01.m
    % Volumetric soil water content at saturation (porosity)
    % (Clapp and Hornberger. 1978. Water Resources Research 14:601-604)
    
    ground.STATVAR.vegetation.soilvar.watsat(i) = 0.8; %, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482];
    %     ground.STATVAR.vegetation.soilvar.soil_texture = 5;       % Soil texture class
    
    h2osoi_sat = ground.STATVAR.vegetation.soilvar.watsat(i) * ground.STATVAR.vegetation.physcon.rhowat * ground.STATVAR.vegetation.soilvar.dz(i);
    
    
    % ground.STATVAR.vegetation.soilvar.soil_water_matric_potential = -50; %should be calculated from soil water content
    
    % ground.STATVAR.vegetation.soilvar.soil_water_matric_potential = -(1/alpha)*(((0w - 0wr) / (0ws - 0wr))^(n /(1-n))-1)^(1/n);
    
    %soil parameters sandy clay loam (p.120, bonan book)
    ground.STATVAR.vegetation.soilvar.alpha = 5.9; %cm -> m  0.059;
    ground.STATVAR.vegetation.soilvar.n = 1.48;
    ground.STATVAR.vegetation.soilvar.Owr = ground.STATVAR.vegetation.soilvar.watsat;
    ground.STATVAR.vegetation.soilvar.Ow = ground.STATVAR.vegetation.soilvar.h2osoi_vol;
    ground.STATVAR.vegetation.soilvar.Ows = 0;
    
    ground.STATVAR.vegetation.soilvar.soil_water_matric_potential = -(1/ground.STATVAR.vegetation.soilvar.alpha)*(((ground.STATVAR.vegetation.soilvar.Ow - ground.STATVAR.vegetation.soilvar.Owr) / (ground.STATVAR.vegetation.soilvar.Ows - ground.STATVAR.vegetation.soilvar.Owr))^(ground.STATVAR.vegetation.soilvar.n /(1-ground.STATVAR.vegetation.soilvar.n)-1))^(1/ground.STATVAR.vegetation.soilvar.n); %m
    
    % Actual water content is some fraction of saturation. These are only used for soil
    % thermal properties and phase change. Note the inconsistency with the use of soil
    % water in the bucket model to calculate the soil wetness factor.
    
    satfrac = 0.2;
    if (ground.STATVAR.vegetation.soilvar.tsoi(i) > ground.STATVAR.vegetation.physcon.tfrz)
        ground.STATVAR.vegetation.soilvar.h2osoi_ice(i) = 0;
%         ground.STATVAR.vegetation.soilvar.h2osoi_liq(i) = satfrac * h2osoi_sat;
    else
%         ground.STATVAR.vegetation.soilvar.h2osoi_liq(i) = 0;
        ground.STATVAR.vegetation.soilvar.h2osoi_ice(i) = satfrac * h2osoi_sat;
    end
    
end
end

