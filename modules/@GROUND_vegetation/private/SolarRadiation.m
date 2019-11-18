
function [vegetation] = SolarRadiation (vegetation)

% %DESCRIPTION:
% Solar radiation transfer through canopy

% %ARGUMENTS:

%    %LOCAL VARIABLES:
%     integer  :: f                                        % Filter index
%     integer  :: p                                        % Patch index for CLM g/l/c/p hierarchy
%     integer  :: c                                        % Column index for CLM g/l/c/p hierarchy
%     integer  :: ic                                       % Aboveground layer index
%     integer  :: ib                                       % Waveband index
%     integer  :: lev                                      % Canopy level index
%     integer  :: j                                        % Sky angle index
%     real(r8) :: angle                                    % Sky angles (5, 15, 25, 35, 45, 55, 65, 75, 85 degrees)
%     real(r8) :: gdirj                                    % Relative projected area of leaf elements in the direction of sky angle
%     real(r8) :: cumpai                                   % Cumulative plant area index (m2/m2)
%     real(r8) :: chil(bounds%begp:bounds%endp)            % Departure of leaf angle from spherical orientation (-0.4 <= xl <= 0.6)
%     real(r8) :: phi1(bounds%begp:bounds%endp)            % Term in Ross-Goudriaan function for gdir
%     real(r8) :: phi2(bounds%begp:bounds%endp)            % Term in Ross-Goudriaan function for gdir
%     real(r8) :: gdir(bounds%begp:bounds%endp)            % Relative projected area of leaf elements in the direction of solar beam
%     real(r8) :: kb(bounds%begp:bounds%endp)              % Direct beam extinction coefficient
%     real(r8) :: wl(bounds%begp:bounds%endp)              % Leaf fraction of canopy
%     real(r8) :: ws(bounds%begp:bounds%endp)              % Stem fraction of canopy
%     real(r8) :: rho(bounds%begp:bounds%endp,1:numrad)    % Leaf/stem reflectance
%     real(r8) :: tau(bounds%begp:bounds%endp,1:numrad)    % Leaf/stem transmittance
%     real(r8) :: omega(bounds%begp:bounds%endp,1:numrad)  % Leaf/stem scattering coefficient
%
%     % For Norman radiation
%     real(r8) :: tbj(bounds%begp:bounds%endp,0:nlevcan)   % Exponential transmittance of direct beam onto canopy layer
%     real(r8) :: tb(bounds%begp:bounds%endp,1:nlevcan)    % Exponential transmittance of direct beam through a single leaf layer

% *** Input ***
xl         = 0.01; % 0.01;                  % Departure of leaf angle from spherical orientation (-)
clump_fac  = 0.7; %1. ;                     % https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010GB003996  % pft, 1 in matlab exercise Bonan (sp_14_03.m)        ; % Foliage clumping index (-)
rhol       = vegetation.leaf.rhol              ; % Leaf reflectance (-)
taul       = vegetation.leaf.taul              ; % Leaf transmittance (-)
rhos       = vegetation.leaf.rhos              ; % Stem reflectance (-)
taus       = vegetation.leaf.taus              ; % Stem transmittance (-)

% vegetation.flux.albsoib = [0.8,0.8]; % Direct beam albedo of ground (soil)
% vegetation.flux.albsoid = [0.8,0.8]; % Diffuse albedo of ground (soil)

% lai        = vegetation.canopy.lai        ; % Leaf area index of canopy (m2/m2)
% sai        = vegetation.mlcanopyinst.sai        ; % Stem area index of canopy (m2/m2)
% ncan       = vegetation.mlcanopyinst.ncan       ; % Number of aboveground layers
% nbot       = vegetation.canopy.nbot      ; % Index for bottom leaf layer
% ntop       = vegetation.canopy.ntop       ; % Index for top leaf layer
% dpai       = vegetation.canopy.dpai       ; % Layer plant area index (m2/m2)
% % sumpai     = vegetation.mlcanopyinst.sumpai     ; % Cumulative plant area index (m2/m2)
% solar_zen  = vegetation.atmos.solar_zenith  ; % Solar zenith angle (radians)

% *** Output ***
% fracsun    = vegetation.flux.fracsun    ; % Sunlit fraction of canopy layer
% fracsha    = vegetation.flux.fracsha    ; % Shaded fraction of canopy layer
% swleaf     = vegetation.flux.swleaf     ; % Leaf absorbed solar radiation (W/m2 leaf)
% apar       = vegetation.flux.apar       ; % Leaf absorbed PAR (umol photon/m2 leaf/s)
% td         = vegetation.flux.td         ;  % Exponential transmittance of diffuse radiation through a single leaf layer

%---------------------------------------------------------------------
% Weight reflectance and transmittance by lai and sai and calculate
% leaf scattering coefficient
%---------------------------------------------------------------------

for f = 1:vegetation.canopy.num_exposedvegp
    p = vegetation.canopy.filter_exposedvegp(f);
    if ((vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p)) > 0.)
        wl(p) = vegetation.canopy.lai(p) / (vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p));
        ws(p) = vegetation.mlcanopyinst.sai(p) / (vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p));
    else
        wl(p) = 0. ;
        ws(p) = 0. ;
    end
end

for ib = 1:vegetation.mlcanopyinst.numrad
    for f = 1:vegetation.canopy.num_exposedvegp
        p = vegetation.canopy.filter_exposedvegp(f);
        rho(p,ib) = max(rhol(p,ib)*wl(p) + rhos(p,ib)*ws(p), 1.e-06 );
        tau(p,ib) = max(taul(p,ib)*wl(p) + taus(p,ib)*ws(p), 1.e-06 );
        omega(p,ib) = rho(p,ib) + tau(p,ib);
    end
end

%---------------------------------------------------------------------
% Direct beam extinction coefficient
%---------------------------------------------------------------------

for f = 1:vegetation.canopy.num_exposedvegp
    p = vegetation.canopy.filter_exposedvegp(f);
    
    chil(p) = min(max(xl(p), -0.4), 0.6);
    if (abs(chil(p)) <= 0.01)
        chil(p) = 0.01 ;
    end
    phi1(p) = 0.5  - 0.633 *chil(p) - 0.330 *chil(p)*chil(p);
    phi2(p) = 0.877  * (1.  - 2. *phi1(p));
    
%     phi1(p) = 0. ;       % Horizontal leaves
%     phi2(p) = 1. ;       % Horizontal leaves
    
    gdir(p) = phi1(p) + phi2(p) * cos(vegetation.sun.solar_zenit(p));
    kb(p) = gdir(p) / cos(vegetation.sun.solar_zenit(p));
    kb(p) = min(kb(p), 40. );
end

%---------------------------------------------------------------------
% Sunlit and shaded fraction of leaf layer
%---------------------------------------------------------------------

for f = 1:vegetation.canopy.num_exposedvegp
    p = vegetation.canopy.filter_exposedvegp(f);
    
    % Zero out for all layers
    
    for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran 1:ncan
        vegetation.flux.fracsun(p,ic) = 0. ;
        vegetation.flux.fracsha(p,ic) = 0. ;
    end
    
    % Leaf layers
    
    for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
        vegetation.flux.fracsun(p,ic) = clump_fac(p) * exp(-kb(p) * vegetation.mlcanopyinst.sumpai(p,ic) * clump_fac(p));
        vegetation.flux.fracsha(p,ic) = 1.  - vegetation.flux.fracsun(p,ic);
    end
    
end

%---------------------------------------------------------------------
% Diffuse transmittance for a single layer estimated for nine sky
% angles in increments of 10 degrees (also needed for longwave radiation)
%---------------------------------------------------------------------

for f = 1:vegetation.canopy.num_exposedvegp
    p = vegetation.canopy.filter_exposedvegp(f);
    
    % Zero out for all layers
    
    for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran 1:ncan
        td(p,ic) = 0.;
    end
    
    % Leaf layers
    
    for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
        for j = 1:9
            angle = (5.  + (j - 1) * 10. ) * pi / 180. ;
            gdirj = phi1(p) + phi2(p) * cos(angle);
            td(p,ic) = td(p,ic) + exp(-gdirj / cos(angle) * vegetation.canopy.dpai(p,ic) * clump_fac(p)) * sin(angle) * cos(angle);
        end
        td(p,ic) = td(p,ic) * 2.  * (10.  * pi / 180. );
    end
    
end

%---------------------------------------------------------------------
% Parameters for Norman radiative transfer
%---------------------------------------------------------------------

% Direct beam transmittance (tb) for a single leaf layer

for f = 1:vegetation.canopy.num_exposedvegp
    p = vegetation.canopy.filter_exposedvegp(f);
    
    % Zero out for all layers
    
    for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran 1:ncan
        tb(p,ic) = 0. ;
    end
    
    % Leaf layers
    
    for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
        tb(p,ic) = exp(-kb(p) * vegetation.canopy.dpai(p,ic) * clump_fac(p));
    end
    
end

% Direct beam transmittance (tbj) uses cumulative plant area index
% above layer j to give unscattered direct beam onto layer j

for f = 1:vegetation.canopy.num_exposedvegp
    p = vegetation.canopy.filter_exposedvegp(f);
    
    % Zero out for all layers
    
    for ic = 1:vegetation.mlcanopyinst.ncan(p) %ic = 0:ncan(p)
        tbj(p,ic) = 0.;
    end
    
    % Leaf layers and ground
    
    cumpai = 0. ;
    tbj(p,vegetation.canopy.ntop(p)) = 1.;
    for ic = vegetation.canopy.ntop(p):-1:vegetation.canopy.nbot(p)
        cumpai = cumpai + vegetation.canopy.dpai(p,ic);
        if (ic > vegetation.canopy.nbot(p))
            lev = ic - 1 ;  % transmittance onto leaf layer below
        elseif (ic == vegetation.canopy.nbot(p))
            lev = 1     ;   % transmittance onto ground % lev = 0
        end
        tbj(p,lev) = exp(-kb(p) * cumpai * clump_fac(p));
    end

end
% % % ---------------------------------------------------------------------
% % % Parameters for Goudriaan radiative transfer
% % % ---------------------------------------------------------------------
% % % 
% % % Diffuse extinction coefficient for canopy
% % % 
% % % for f = 1:vegetation.canopy.num_exposedvegp
% % %     p = vegetation.canopy.filter_exposedvegp(f);
% % %     kd(p) = 0. ;
% % %     for j = 1:9
% % %         angle = (5.  + (j - 1) * 10. ) * pi / 180. ;
% % %         gdirj = phi1(p) + phi2(p) * cos(angle);
% % %         kd(p) = kd(p) + exp(-gdirj / cos(angle) * (vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p))) * sin(angle) * cos(angle);
% % %     end
% % %     kd(p) = kd(p) * 2.  * (10.  * pi / 180. );
% % %     if ((vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p)) > 0. )
% % %         kd(p) = -log(kd(p)) / (vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p));
% % %     else
% % %         kd(p) = 0. ;
% % %     end
% % % end
% % % 
% % % Adjust extinction coefficents for scattering and calculate canopy albefors
% % % 
% % % for ib = 1:vegetation.mlcanopyinst.numrad
% % %     for f = 1:vegetation.canopy.num_exposedvegp
% % %         p = vegetation.canopy.filter_exposedvegp(f);
% % %         c = p;
% % %         
% % %         Adjust for scattering
% % %         
% % %         kbm(p,ib) = kb(p) * sqrt(1.  - omega(p,ib));
% % %         kdm(p,ib) = kd(p) * sqrt(1.  - omega(p,ib));
% % %         
% % %         Vegetation albefor, horizontal leaves
% % %         
% % %         albvegh = (1.  - sqrt(1.  - omega(p,ib))) / (1.  + sqrt(1.  - omega(p,ib)));
% % %         
% % %         Direct beam vegetation albefor, non-horizontal leaves
% % %         
% % %         albvegb = 2.  * kb(p) / (kb(p) + kd(p)) * albvegh;
% % %         
% % %         Diffuse vegetation albefor, non-horizontal leaves
% % %         
% % %         albvegd = 0.;
% % %         for j = 1:9
% % %             angle = (5.  + (j - 1) * 10. ) * pi / 180. ;
% % %             gdirj = phi1(p) + phi2(p) * cos(angle);
% % %             kbj = gdirj / cos(angle);
% % %             albvegbj = 2.  * kbj / (kbj + kd(p)) * albvegh;
% % %             albvegd = albvegd + albvegbj * sin(angle) * cos(angle);
% % %         end
% % %         albvegd = albvegd * 2.  * (10.  * pi / 180. );
% % %         
% % %         Effective canopy albefor, including soil
% % %         
% % %         albcanb(p,ib) = albvegb + (vegetation.flux.albsoib(c,ib) - albvegb) ...
% % %             * exp(-2.  * kbm(p,ib) * (vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p)) * clump_fac(p));
% % %         albcand(p,ib) = albvegd + (vegetation.flux.albsoid(c,ib) - albvegd) ...
% % %             * exp(-2.  * kdm(p,ib) * (vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p)) * clump_fac(p));
% % %         
% % %     end
% % % end
% % % 

%---------------------------------------------------------------------
% Parameters for two-stream radiative transfer
%---------------------------------------------------------------------

% avmu - average inverse diffuse optical depth per unit leaf area
% 
% for f = 1:vegetation.canopy.num_exposedvegp
%     p = vegetation.canopy.filter_exposedvegp(f);
%     avmu(p) = ( 1.  - phi1(p)/phi2(p) * log((phi1(p)+phi2(p))/phi1(p)) ) / phi2(p);
% end
% 
% % betad - upscatter parameter for diffuse radiation
% % betab - upscatter parameter for direct beam radiation
% 
% for ib = 1:vegetation.mlcanopyinst.numrad
%     for f = 1:vegetation.canopy.num_exposedvegp
%         p = vegetation.canopy.filter_exposedvegp(f);
%         
%         % upscatter parameter for diffuse radiation
%         
%         betad(p,ib) = 0.5  / omega(p,ib)...
%             * (rho(p,ib) + tau(p,ib) + (rho(p,ib)-tau(p,ib)) * ((1. +chil(p))/2. )^2 );
%         
%         % upscatter parameter for direct beam radiation
%         
%         tmp0 = gdir(p) + phi2(p) * cos(vegetation.sun.solar_elevation(p));
%         tmp1 = phi1(p) * cos(vegetation.sun.solar_elevation(p));
%         tmp2 = 1.  - tmp1/tmp0 * log((tmp1+tmp0)/tmp1);
%         asu = 0.5  * omega(p,ib) * gdir(p) / tmp0 * tmp2;
%         betab(p,ib) = (1.  + avmu(p)*kb(p)) / (omega(p,ib)*avmu(p)*kb(p)) * asu;
%         
%     end
% end

%---------------------------------------------------------------------
% Calculate radiative transfer through canopy
%---------------------------------------------------------------------
[vegetation] = NormanRadiation (rho, tau, omega, td, tb, tbj, vegetation); %, p, ic
% [vegetation] = GoudriaanRadiation (omega, kb, kbm, kdm, albcanb, albcand, vegetation);
% [vegetation] = TwoStreamRadiation (omega, kb, avmu, betad, betab, vegetation);


%     if (light == 1)
%        call NormanRadiation (bounds, vegetation.physcon.num_exposedvegp, vegetation.physcon.filter_exposedvegp, &
%        rho, tau, omega, tb, td, tbj, surfalb_inst, mlcanopy_inst)
%     else if (light == 2)
%        call GoudriaanRadiation (bounds, vegetation.physcon.num_exposedvegp, vegetation.physcon.filter_exposedvegp, &
%        omega, kb, kbm, kdm, albcanb, albcand, mlcanopy_inst)
%     else if (light == 3)
%        call TwoStreamRadiation (bounds, vegetation.physcon.num_exposedvegp, vegetation.physcon.filter_exposedvegp, &
%        omega, kb, avmu, betad, betab, surfalb_inst, mlcanopy_inst)
%     end

% APAR per unit sunlit and shaded leaf area

for f = 1:vegetation.canopy.num_exposedvegp
    p = vegetation.canopy.filter_exposedvegp(f);
    for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran 1:ncan
        vegetation.flux.apar(p,ic,vegetation.params.sun) = vegetation.flux.swleaf(p,ic,vegetation.params.sun,vegetation.params.vis) * 4.6 ;
        vegetation.flux.apar(p,ic,vegetation.params.sha) = vegetation.flux.swleaf(p,ic,vegetation.params.sha,vegetation.params.vis) * 4.6 ;
    end
end
% *** Output ***
% vegetation.flux.fracsun    = fracsun    ; % Sunlit fraction of canopy layer
% vegetation.flux.fracsha    = fracsha    ; % Shaded fraction of canopy layer
% vegetation.flux.swleaf     = swleaf     ; % Leaf absorbed solar radiation (W/m2 leaf)
% vegetation.flux.apar       = apar       ; % Leaf absorbed PAR (umol photon/m2 leaf/s)
vegetation.flux.td = td         ;  % Exponential transmittance of diffuse radiation through a single leaf layer

end
