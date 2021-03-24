function [vegetation] = SoilResistance(vegetation)
%(num_exposedvegp, filter_exposedvegp, soilstate_inst, waterstate_inst, mlcanopy_inst, clm_varpar, clm_varcon, pftcon, ColumnType, PatchType, SoilStateType, WaterStateType, CanopyFluxesMultilayerType)
%
% %DESCRIPTION:
% Calculate soil hydraulic resistance and water uptake from each soil layer

%LOCAL VARIABLES:
%      f                                % Filter index
%      p                                % Patch index for CLM g/l/c/p hierarchy
%      c                                % Column index for CLM g/l/c/p hierarchy
%      j                                % Soil layer index
%      s                                % Soil layer water content relative to saturation (fraction)
%      root_cross_sec_area              % Root cross-sectional area (m2 root)
%      root_biomass_density             % Root biomass density (g biomass / m3 soil)
%      root_length_density              % Root length density (m root / m3 soil)
%      root_dist                        % Mean distance between roots (m)
%      hk                               % Hydraulic conductivity (mm/s -> mmol/m/s/MPa)
%      soilr1                           % Soil-to-root resistance (MPa.s.m2/mmol H2O)
%      soilr2                           % Root-to-stem resistance (MPa.s.m2/mmol H2O)
%      soilr                            % Belowground resistance (MPa.s.m2/mmol H2O)
%      smp_mpa(nlevsoi)                 % Soil matric potential (MPa)
%      evap(nlevsoi)                    % Maximum transpiration (mmol H2O/m2/s)
%      totevap                          % Total maximum transpiration (mmol H2O/m2/s)

vegetation.physcon.head = vegetation.physcon.denh20*vegetation.physcon.grav*1.e-06;
%---------------------------------------------------------------------
% minlwp       = vegetation.pftcon.minlwp                 ;  % Minimum leaf water potential (MPa)
% root_radius  = vegetation.pftcon.root_radius            ;  % Fine root radius (m)
% root_density = vegetation.pftcon.root_density           ;  % Fine root density (g biomass / m3 root)
% root_resist  = vegetation.pftcon.root_resist            ;  % Hydraulic resistivity of root tissue (MPa.s.g/mmol H2O)

% values from here: https://slideplayer.com/slide/12736113/76/images/21/Clapp+and+Hornberger+(1978)+parameters+for+Soil+Moisture+Characteristic+functions+based+on+analysis+of+1845+soils.+Values+in+parentheses+are+standard+deviations..jpg
% dz           = 1; %vegetation.soilvar.dz         ;  % Soil layer thickness (m)
% watsat       = 0.4; %vegetation.soilvar.watsat     ;  % Soil layer volumetric water content at saturation (porosity)
% hksat        = 3.; %soilstate_inst.hksat_col      ;  % Soil layer hydraulic conductivity at saturation (mm H2O/s)
% bsw          = 4; %soilstate_inst.bsw_col        ;  % Soil layer Clapp and Hornberger "b" parameter
% smp_l        = 1; %soilstate_inst.smp_l_col      ;  % Soil layer matric potential (mm)
% rootfr       = 0.5; %soilstate_inst.rootfr_patch   ;  % Fraction of roots in each soil layer
% h2osoi_vol   = 0.1; %waterstate_inst.h2osoi_vol_col;  % Soil layer volumetric water content (m3/m3)
% h2osoi_ice   = 0; %waterstate_inst.h2osoi_ice_col;  % Soil layer ice lens (kg/m2)

% Soil and root resistances for each layer

for f = 1:vegetation.canopy.num_exposedvegp
    p = vegetation.canopy.filter_exposedvegp(f);
    c = p;
    root_cross_sec_area = pi .* vegetation.pftcon.root_radius(p)^2;
    
    vegetation.mlcanopyinst.rsoil(p) = 0.;
    for j = 1:vegetation.soilvar.nsoi
        
        %get variables from GROUND module 
        hk = vegetation.PARENT_GROUND.STATVAR.hydraulicConductivity(j,1) .*1000;
        smp_mpa(j) = vegetation.PARENT_GROUND.STATVAR.waterPotential(j,1) .* vegetation.physcon.head;
        vegetation.soilvar.dz(c,j) = vegetation.PARENT_GROUND.STATVAR.layerThick(j,1);
        vegetation.soilvar.h2osoi_ice(c,j) = vegetation.PARENT_GROUND.STATVAR.ice(j,1) ./ vegetation.PARENT_GROUND.STATVAR.layerThick(j,1) ./ vegetation.PARENT_GROUND.STATVAR.area(j,1);
        
        %SEBAS:remove + add
       %  s = max(min(vegetation.soilvar.h2osoi_vol(c,j)/vegetation.soilvar.watsat(c,j), 1.), 0.01);
       % hk = vegetation.soilvar.hksat(c,j) .* s^(2  .* vegetation.soilvar.bsw(c,j) + 3.);           % mm/s
        %hk = vegetation.soilvar.hk(c,j) .* 1000;  %read from GROUND in m/sec
        %end SEBAS:remove + add
        
        hk = hk .* 1.e-03  ./ vegetation.physcon.head;                                % mm/s -> m/s -> m2/s/MPa
        hk = hk .* vegetation.physcon.denh20 ./ vegetation.physcon.mmh2o .* 1000.;                        % m2/s/MPa -> mmol/m/s/MPa

        %für smp_l van genuchten formel benutzen n, m, parameter p.119
        
%         vegetation.soilvar.smp_l = vegetation.soilvar.soil_water_matric_potential*1000.;
%         smp_mpa(j) = vegetation.soilvar.smp_l(c,j) .* 1.e-03  .* vegetation.physcon.head;                % mm -> m -> MPa
        
        % Root biomass density: g biomass / m3 soil
                
        root_biomass_density = vegetation.mlcanopyinst.root_biomass(p) .* vegetation.soilvar.rootfr(p,j) ./ vegetation.soilvar.dz(c,j);
        root_biomass_density = max(root_biomass_density, 1.e-10 );
        
        
        % Root length density: m root per m3 soil
        
        root_length_density = root_biomass_density ./ (vegetation.pftcon.root_density(p) .* root_cross_sec_area);
        
        % Distance between roots: m
        
        root_dist = sqrt (1. ./ (root_length_density .* pi));
        
        % Soil-to-root resistance (MPa.s.m2/mmol H2O)
        
%         root_length_density
        
        soilr1 = log(root_dist./vegetation.pftcon.root_radius(p)) ./ (2.  .* pi .* root_length_density .* vegetation.soilvar.dz(c,j) .* hk);
        
        % Root-to-stem resistance (MPa.s.m2/mmol H2O)
        
        soilr2 = vegetation.pftcon.root_resist(p) ./ (root_biomass_density .* vegetation.soilvar.dz(c,j));
        
        % Belowground resistance (MPa.s.m2/mmol H2O)
        
         %disp(soilr1)
        
        soilr = soilr1 + soilr2;
        
        % Total belowground resistance. First sum the conductances (1/soilr)
        % for each soil layer and then convert back to a resistance after the
        % summation.
        
        %SEBAS CHANGE: set to zero of there is ice in a layer, just as below!!
        %vegetation.mlcanopyinst.rsoil(p) = vegetation.mlcanopyinst.rsoil(p) + 1. ./ soilr;
        vegetation.mlcanopyinst.rsoil(p) = vegetation.mlcanopyinst.rsoil(p) + double(vegetation.soilvar.h2osoi_ice(c,j) == 0) ./ soilr;
        %END CHANGE
        
        % Maximum transpiration for each layer (mmol H2O/m2/s). No negative
        % transpiration and no transpiration from frozen soil.
        
        evap(j) = (smp_mpa(j) - vegetation.leaf.minlwp(p)) ./ soilr;
        evap(j) = max (evap(j), 0.);
        if (vegetation.soilvar.h2osoi_ice(c,j) > 0.)
            evap(j) = 0 ;
        end
        
        vegetation.soilvar.transp_per_layer(j) = evap(j)/1000; % Simone: mmol -> mol
    end
    
    % Belowground resistance: resistance = 1 / conductance
    
    vegetation.mlcanopyinst.rsoil(p) = vegetation.canopy.lai(p) ./ vegetation.mlcanopyinst.rsoil(p);
    %SEBAST CHANGE: set to some high value of entire soil is frozen
    if vegetation.mlcanopyinst.rsoil(p) == Inf
        vegetation.mlcanopyinst.rsoil(p) = 1e20;
    end
    % END CHNAGE
    
    % Weighted soil water potential (MPa) and fractional uptake from soil layers
    
    totevap = sum(evap);
    vegetation.mlcanopyinst.psis(p) = 0;
    vegetation.mlcanopyinst.soil_et_loss(p,:) = 0;
    
    for j = 1:vegetation.soilvar.nsoi
        vegetation.mlcanopyinst.psis(p) = vegetation.mlcanopyinst.psis(p) + smp_mpa(j) .* evap(j);
        if (totevap > 0.)
            vegetation.mlcanopyinst.soil_et_loss(p,j) = evap(j) ./ totevap;
        else
            %SEBAS CHANGED
            %vegetation.mlcanopyinst.soil_et_loss(p,j) = 1. ./ vegetation.mlcanopyinst.nlevgrnd;  %SEBAS: THIS IS 1, WHY?
            vegetation.mlcanopyinst.soil_et_loss(p,j) = 0; %1. ./ vegetation.mlcanopyinst.nlevgrnd;  %SEBAS: THIS IS 1, WHY?
            %END CHANGE
        end
    end
    
    if (totevap > 0.)
        vegetation.mlcanopyinst.psis(p) = vegetation.mlcanopyinst.psis(p) ./ totevap;
    else
        vegetation.mlcanopyinst.psis(p) = vegetation.leaf.minlwp(p);
    end
end
end
