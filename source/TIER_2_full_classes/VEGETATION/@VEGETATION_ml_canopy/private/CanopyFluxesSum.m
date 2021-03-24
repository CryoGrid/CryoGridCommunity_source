
function [vegetation] = CanopyFluxesSum (vegetation, p) %, ic, il
%
% %DESCRIPTION:
% Sum leaf and soil fluxes
%
% %USES:
%     turb_type
ivis = vegetation.params.vis; % Array index for visible waveband
inir = vegetation.params.nir; % Array index for near-infrared waveband
isun = vegetation.params.sun; % Array index for sunlit leaf
isha = vegetation.params.sha; % Array index for shaded leaf

% %LOCAL VARIABLES:
%     ic                  % Aboveground layer index
%     err                 % Energy imbalance
%     radin               % Incoming radiation
%     raforut             % Outgoing radiation
%     avail               % Available energy
%     flux                % Turbulent fluxes + storage
%     fracgreen           % Green fraction of plant area index: lai/(lai+sai)
%---------------------------------------------------------------------

% vegetation.mlcanopyinst.ncan        = vegetation.mlcanopyinst.ncan         ; % Number of aboveground layers
% vegetation.mlcanopyinst.nbot        = vegetation.mlcanopyinst.nbot         ; % Index for bottom leaf layer
% vegetation.canopy.ntop        = vegetation.canopy.ntop         ; % Index for top leaf layer
% vegetation.mlcanopyinst.dpai        = vegetation.mlcanopyinst.dpai         ; % Layer plant area index (m2/m2)
% vegetation.mlcanopyinst.tref        = vegetation.mlcanopyinst.tref         ; % Air temperature at reference height (K)
% vegetation.mlcanopyinst.fwet        = vegetation.mlcanopyinst.fwet         ; % Fraction of plant area index that is wet
% vegetation.mlcanopyinst.fdry        = vegetation.mlcanopyinst.fdry         ; % Fraction of plant area index that is green and dry
% vegetation.mlcanopyinst.shair       = vegetation.mlcanopyinst.shair        ; % Canopy air sensible heat flux (W/m2)
% vegetation.mlcanopyinst.etair       = vegetation.mlcanopyinst.etair        ; % Canopy air water vapor flux (mol H2O/m2/s)
% vegetation.mlcanopyinst.stair       = vegetation.mlcanopyinst.stair        ; % Canopy air storage heat flux (W/m2)
% vegetation.mlcanopyinst.irleaf      = vegetation.mlcanopyinst.irleaf       ; % Leaf absorbed longwave radiation (W/m2 leaf)
% vegetation.mlcanopyinst.rnleaf      = vegetation.mlcanopyinst.rnleaf       ; % Leaf net radiation (W/m2 leaf)
% vegetation.mlcanopyinst.stleaf      = vegetation.mlcanopyinst.stleaf       ; % Leaf storage heat flux (W/m2 leaf)
% vegetation.mlcanopyinst.shleaf      = vegetation.mlcanopyinst.shleaf       ; % Leaf sensible heat flux (W/m2 leaf)
% vegetation.mlcanopyinst.lhleaf      = vegetation.mlcanopyinst.lhleaf       ; % Leaf latent heat flux (W/m2 leaf)
% vegetation.mlcanopyinst.trleaf      = vegetation.mlcanopyinst.trleaf       ; % Leaf transpiration flux (mol H2O/m2 leaf/s)
% vegetation.mlcanopyinst.evleaf      = vegetation.mlcanopyinst.evleaf       ; % Leaf evaporation flux (mol H2O/m2 leaf/s)
% vegetation.flux.swleaf              = vegetation.flux.swleaf       ; % Leaf absorbed solar radiation (W/m2 leaf)
% vegetation.mlcanopyinst.an          = vegetation.mlcanopyinst.an           ; % Leaf net photosynthesis (umol CO2/m2 leaf/s)
% vegetation.mlcanopyinst.ag          = vegetation.mlcanopyinst.ag           ; % Leaf gross photosynthesis (umol CO2/m2 leaf/s)
% vegetation.flux.fracsun             = vegetation.flux.fracsun      ; % Sunlit fraction of canopy layer
% vegetation.flux.fracsha             = vegetation.flux.fracsha      ; % Shaded fraction of canopy layer
% vegetation.flux.swveg               = vegetation.flux.swveg        ; % Absorbed solar radiation, vegetation (W/m2)
% vegetation.mlcanopyinst.rnsoi       = vegetation.mlcanopyinst.rnsoi        ; % Net radiation, ground (W/m2)
% vegetation.mlcanopyinst.shsoi       = vegetation.mlcanopyinst.shsoi        ; % Sensible heat flux, ground (W/m2)
% vegetation.mlcanopyinst.lhsoi       = vegetation.mlcanopyinst.lhsoi        ; % Latent heat flux, ground (W/m2)
% vegetation.mlcanopyinst.gsoi        = vegetation.mlcanopyinst.gsoi         ; % Soil heat flux (W/m2)
% vegetation.flux.swsoi               = vegetation.flux.swsoi        ; % Absorbed solar radiation, ground (W/m2)
% vegetation.flux.irsoi               = vegetation.flux.irsoi        ; % Absorbed longwave radiation, ground (W/m2)
% vegetation.mlcanopyinst.etsoi       = vegetation.mlcanopyinst.etsoi        ; % Water vapor flux, ground (mol H2O/m2/s)
% vegetation.atmos.swskyb             = vegetation.atmos.swskyb       ; % Atmospheric direct beam solar radiation (W/m2)
% vegetation.atmos.swskyd             = vegetation.atmos.swskyd       ; % Atmospheric diffuse solar radiation (W/m2)
% vegetation.mlcanopyinst.irsky       = vegetation.mlcanopyinst.irsky        ; % Atmospheric longwave radiation (W/m2)
% vegetation.flux.albcan              = vegetation.flux.albcan       ; % Albefor above canopy
% vegetation.mlcanopyinst.ircan       = vegetation.mlcanopyinst.ircan        ; % Upward longwave radiation above canopy (W/m2)
% *** Output ***
% vegetation.mlcanopyinst.sw_prof     = vegetation.mlcanopyinst.sw_prof      ; % Canopy layer absorbed solar radiation (W/m2)
% vegetation.mlcanopyinst.ir_prof     = vegetation.mlcanopyinst.ir_prof      ; % Canopy layer absorbed longwave radiation (W/m2)
% vegetation.mlcanopyinst.rn_prof     = vegetation.mlcanopyinst.rn_prof      ; % Canopy layer net radiation (W/m2)
% vegetation.mlcanopyinst.st_prof     = vegetation.mlcanopyinst.st_prof      ; % Canopy layer storage heat flux (W/m2)
% vegetation.mlcanopyinst.sh_prof     = vegetation.mlcanopyinst.sh_prof      ; % Canopy layer sensible heat flux (W/m2)
% vegetation.mlcanopyinst.sh_prof     = vegetation.mlcanopyinst.lh_prof      ; % Canopy layer latent heat flux (W/m2)
% vegetation.mlcanopyinst.et_prof     = vegetation.mlcanopyinst.et_prof      ; % Canopy layer water vapor flux (mol H2O/m2/s)
% vegetation.mlcanopyinst.fc_prof     = vegetation.mlcanopyinst.fc_prof      ; % Canopy layer CO2 flux (umol CO2/m2/s)
% vegetation.mlcanopyinst.rnet        = vegetation.mlcanopyinst.rnet         ; % Net radiation (W/m2)
% vegetation.mlcanopyinst.stflx       = vegetation.mlcanopyinst.stflx        ; % Canopy storage heat flux (W/m2)
% vegetation.mlcanopyinst.shflx       = vegetation.mlcanopyinst.shflx        ; % Sensible heat flux (W/m2)
% vegetation.mlcanopyinst.lhflx       = vegetation.mlcanopyinst.lhflx        ; % Latent heat flux (W/m2)
% vegetation.mlcanopyinst.etflx       = vegetation.mlcanopyinst.etflx        ; % Water vapor flux (mol H2O/m2/s)
% vegetation.flux.irveg               = vegetation.flux.irveg        ; % Absorbed longwave radiation, vegetation (W/m2)
% vegetation.mlcanopyinst.irvegsun    = vegetation.mlcanopyinst.irvegsun     ; % Absorbed longwave radiation, sunlit canopy (W/m2)
% vegetation.mlcanopyinst.irvegsha    = vegetation.mlcanopyinst.irvegsha     ; % Absorbed longwave radiation, shaded canopy (W/m2)
% vegetation.mlcanopyinst.shveg       = vegetation.mlcanopyinst.shveg        ; % Sensible heat flux, vegetation (W/m2)
% vegetation.mlcanopyinst.shvegsun    = vegetation.mlcanopyinst.shvegsun     ; % Sensible heat flux, sunlit canopy (W/m2)
% vegetation.mlcanopyinst.shvegsha    = vegetation.mlcanopyinst.shvegsha     ; % Sensible heat flux, shaded canopy (W/m2)
% vegetation.mlcanopyinst.lhveg       = vegetation.mlcanopyinst.lhveg        ; % Latent heat flux, vegetation (W/m2)
% vegetation.mlcanopyinst.lhvegsun    = vegetation.mlcanopyinst.lhvegsun     ; % Latent heat flux, sunlit canopy (W/m2)
% vegetation.mlcanopyinst.lhvegsha    = vegetation.mlcanopyinst.lhvegsha     ; % Latent heat flux, shaded canopy (W/m2)
% vegetation.mlcanopyinst.etveg       = vegetation.mlcanopyinst.etveg        ; % Water vapor flux, vegetation (mol H2O/m2/s)
% vegetation.mlcanopyinst.etvegsun    = vegetation.mlcanopyinst.etvegsun     ; % Water vapor flux, sunlit canopy (mol H2O/m2/s)
% vegetation.mlcanopyinst.etvegsha    = vegetation.mlcanopyinst.etvegsha     ; % Water vapor flux, shaded canopy (mol H2O/m2/s)
% vegetation.mlcanopyinst.gppveg      = vegetation.mlcanopyinst.gppveg       ; % Gross primary production (umol CO2/m2/s)
% vegetation.mlcanopyinst.gppvegsun   = vegetation.mlcanopyinst.gppvegsun    ; % Gross primary production, sunlit canopy (umol CO2/m2/s)
% vegetation.mlcanopyinst.gppvegsha   = vegetation.mlcanopyinst.gppvegsha    ;  % Gross primary production, shaded canopy (umol CO2/m2/s)

% Leaf flux profiles

for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan %(p) %Fortran ic = 1:ncan
    if (vegetation.canopy.dpai(p,ic) > 0.)
        vegetation.mlcanopyinst.ir_prof(p,ic) = vegetation.mlcanopyinst.irleaf(p,ic) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.sw_prof(p,ic,ivis) = (vegetation.flux.swleaf(p,ic,isun,ivis)*vegetation.flux.fracsun(p,ic) + vegetation.flux.swleaf(p,ic,isha,ivis)*vegetation.flux.fracsha(p,ic)) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.sw_prof(p,ic,inir) = (vegetation.flux.swleaf(p,ic,isun,inir)*vegetation.flux.fracsun(p,ic) + vegetation.flux.swleaf(p,ic,isha,inir)*vegetation.flux.fracsha(p,ic)) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.rn_prof(p,ic) = (vegetation.mlcanopyinst.rnleaf(p,ic,isun)*vegetation.flux.fracsun(p,ic) + vegetation.mlcanopyinst.rnleaf(p,ic,isha)*vegetation.flux.fracsha(p,ic)) * vegetation.canopy.dpai(p,ic);
        % canopy layer storage heat flux = leaf storage heat flux sunlit * sunlit canopy fraction + leaf storage heat flux shaded * shaded canopy fraction * layer plant area index
        vegetation.mlcanopyinst.st_prof(p,ic) = (vegetation.mlcanopyinst.stleaf(p,ic,isun)*vegetation.flux.fracsun(p,ic) + vegetation.mlcanopyinst.stleaf(p,ic,isha)*vegetation.flux.fracsha(p,ic)) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.sh_prof(p,ic) = (vegetation.mlcanopyinst.shleaf(p,ic,isun)*vegetation.flux.fracsun(p,ic) + vegetation.mlcanopyinst.shleaf(p,ic,isha)*vegetation.flux.fracsha(p,ic)) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.lh_prof(p,ic) = (vegetation.mlcanopyinst.lhleaf(p,ic,isun)*vegetation.flux.fracsun(p,ic) + vegetation.mlcanopyinst.lhleaf(p,ic,isha)*vegetation.flux.fracsha(p,ic)) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.et_prof(p,ic) = (vegetation.mlcanopyinst.evleaf(p,ic,isun) + vegetation.mlcanopyinst.trleaf(p,ic,isun)) * vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic) ...
            + (vegetation.mlcanopyinst.evleaf(p,ic,isha) + vegetation.mlcanopyinst.trleaf(p,ic,isha)) * vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic);
        fracgreen = vegetation.mlcanopyinst.fdry(p,ic) / (1.  - vegetation.mlcanopyinst.fwet(p,ic));
        vegetation.mlcanopyinst.fc_prof(p,ic) = (vegetation.mlcanopyinst.an(p,ic,isun)*vegetation.flux.fracsun(p,ic) + vegetation.mlcanopyinst.an(p,ic,isha)*vegetation.flux.fracsha(p,ic)) * vegetation.canopy.dpai(p,ic) * fracgreen;
    else
        vegetation.mlcanopyinst.ir_prof(p,ic) = 0.;
        vegetation.mlcanopyinst.sw_prof(p,ic,ivis) = 0.;
        vegetation.mlcanopyinst.sw_prof(p,ic,inir) = 0.;
        vegetation.mlcanopyinst.rn_prof(p,ic) = 0.;
        vegetation.mlcanopyinst.st_prof(p,ic) = 0.;
        vegetation.mlcanopyinst.sh_prof(p,ic) = 0.;
        vegetation.mlcanopyinst.lh_prof(p,ic) = 0.;
        vegetation.mlcanopyinst.et_prof(p,ic) = 0.;
        vegetation.mlcanopyinst.fc_prof(p,ic) = 0.;
    end
end

% Sum leaf fluxes

vegetation.flux.irveg(p) = 0.;
vegetation.mlcanopyinst.stflx(p) = 0.;
vegetation.mlcanopyinst.shveg(p) = 0.;
vegetation.mlcanopyinst.lhveg(p) = 0.;
vegetation.mlcanopyinst.etveg(p) = 0.;
vegetation.mlcanopyinst.gppveg(p) = 0.;

for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
    vegetation.flux.irveg(p) = vegetation.flux.irveg(p) + vegetation.mlcanopyinst.ir_prof(p,ic);
    vegetation.mlcanopyinst.stflx(p) = vegetation.mlcanopyinst.stflx(p) + vegetation.mlcanopyinst.st_prof(p,ic); % st_prof = canopy layer storage heat flux
    vegetation.mlcanopyinst.shveg(p) = vegetation.mlcanopyinst.shveg(p) + vegetation.mlcanopyinst.sh_prof(p,ic); % sh_prof = canopy layer latent heat flux
    vegetation.mlcanopyinst.lhveg(p) = vegetation.mlcanopyinst.lhveg(p) + vegetation.mlcanopyinst.lh_prof(p,ic); % lh_prof = canopy layer latent heat flux
    vegetation.mlcanopyinst.etveg(p) = vegetation.mlcanopyinst.etveg(p) + vegetation.mlcanopyinst.et_prof(p,ic); % et_prof = canopy layer water vapor flux
    if (vegetation.canopy.dpai(p,ic) > 0)
        fracgreen = vegetation.mlcanopyinst.fdry(p,ic) / (1.  - vegetation.mlcanopyinst.fwet(p,ic));
        vegetation.mlcanopyinst.gppveg(p) = vegetation.mlcanopyinst.gppveg(p) + (vegetation.mlcanopyinst.ag(p,ic,isun)*vegetation.flux.fracsun(p,ic) + vegetation.mlcanopyinst.ag(p,ic,isha)*vegetation.flux.fracsha(p,ic)) * vegetation.canopy.dpai(p,ic) * fracgreen;
    end
end

% Check energy balance for conservation

% absorbed solar radiation vegetation vis + absorbed solar radiation vegetation nir + absorbed longwave radiation vegetation - sensible heat flux vegetation - latent heat flux vegetation - canopy storage heat flux
err = vegetation.flux.swveg(p,ivis) + vegetation.flux.swveg(p,inir) + vegetation.flux.irveg(p) - vegetation.mlcanopyinst.shveg(p) - vegetation.mlcanopyinst.lhveg(p) - vegetation.mlcanopyinst.stflx(p);
if (abs(err) >= 1.e-03 )
    disp ('ERROR: CanopyFluxesSum: energy conservation error (1)');
    disp(err);
end

% Soil fluxes
ic = 1; %Fortran ic = 0
vegetation.mlcanopyinst.sw_prof(p,ic,ivis) = vegetation.flux.swsoi(p,ivis);
vegetation.mlcanopyinst.sw_prof(p,ic,inir) = vegetation.flux.swsoi(p,inir);
vegetation.mlcanopyinst.rn_prof(p,ic) = vegetation.mlcanopyinst.rnsoi(p);
vegetation.mlcanopyinst.st_prof(p,ic) = 0;
vegetation.mlcanopyinst.sh_prof(p,ic) = vegetation.mlcanopyinst.shsoi(p);
vegetation.mlcanopyinst.sh_prof(p,ic) = vegetation.mlcanopyinst.lhsoi(p);
vegetation.mlcanopyinst.et_prof(p,ic) = vegetation.mlcanopyinst.etsoi(p);
vegetation.mlcanopyinst.fc_prof(p,ic) = 0;

% Turbulent fluxes

% how is this translated?

% % switch (turb_type)
% % case (0)
% %         % Sum of layer fluxes
% %         vegetation.mlcanopyinst.shflx(p) = vegetation.mlcanopyinst.shveg(p) + shsoi(p);
% %         vegetation.mlcanopyinst.lhflx(p) = vegetation.mlcanopyinst.lhveg(p) + lhsoi(p);
% %         vegetation.mlcanopyinst.etflx(p) = vegetation.mlcanopyinst.etveg(p) + etsoi(p);
% % case (1&2&3&4)
% Turbulent fluxes are at the top of the canopy
vegetation.mlcanopyinst.shflx(p) = vegetation.mlcanopyinst.shair(p,vegetation.canopy.ntop(p));
vegetation.mlcanopyinst.etflx(p) = vegetation.mlcanopyinst.etair(p,vegetation.canopy.ntop(p));
vegetation.mlcanopyinst.lhflx(p) = vegetation.mlcanopyinst.etair(p,vegetation.canopy.ntop(p)) * LatVap(vegetation.mlcanopyinst.tref(p), vegetation);
% %         case default
% %             disp ('ERROR: CanopyFluxesMultilayer: turbulence type not valid');
% % end
% %
% Add canopy air heat storage to storage flux

for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p) %Fortran ic = 1:ncan
    vegetation.mlcanopyinst.stflx(p) = vegetation.mlcanopyinst.stflx(p) + vegetation.mlcanopyinst.stair(p,ic);
end

% Overall energy balance check:
% radiation in - radiation out - soil heat = available energy = turbulent flux + canopy storage flux

vegetation.mlcanopyinst.rnet(p) = vegetation.flux.swveg(p,ivis) + vegetation.flux.swveg(p,inir) + vegetation.flux.swsoi(p,ivis) + vegetation.flux.swsoi(p,inir) + vegetation.flux.irveg(p) + vegetation.flux.irsoi(p);
radin = vegetation.atmos.swskyb(p,ivis) + vegetation.atmos.swskyd(p,ivis) + vegetation.atmos.swskyb(p,inir) + vegetation.atmos.swskyd(p,inir) + vegetation.mlcanopyinst.irsky(p);
radout = vegetation.flux.albcan(p,ivis)*(vegetation.atmos.swskyb(p,ivis)+vegetation.atmos.swskyd(p,ivis)) + vegetation.flux.albcan(p,inir)*(vegetation.atmos.swskyb(p,inir)+vegetation.atmos.swskyd(p,inir)) + vegetation.mlcanopyinst.ircan(p);

err = vegetation.mlcanopyinst.rnet(p) - (radin - radout);
if (abs(err) > 0.01)
    disp ('ERROR: CanopyFluxesSum: energy conservation error (2)');
end

% avail = radin - radout - vegetation.mlcanopyinst.gsoi(p);
% flux = vegetation.mlcanopyinst.shflx(p) + vegetation.mlcanopyinst.lhflx(p) + vegetation.mlcanopyinst.stflx(p);
% err = avail - flux;
% if (abs(err) > 0.01)
%     disp ('ERROR: CanopyFluxesSum: energy conservation error (3)');
%     disp (err);
% end

% Sunlit and shaded canopy fluxes

vegetation.mlcanopyinst.irvegsun(p) = 0.;
vegetation.mlcanopyinst.irvegsha(p) = 0.;
vegetation.mlcanopyinst.shvegsun(p) = 0.;
vegetation.mlcanopyinst.shvegsha(p) = 0.;
vegetation.mlcanopyinst.lhvegsun(p) = 0.;
vegetation.mlcanopyinst.lhvegsha(p) = 0.;
vegetation.mlcanopyinst.etvegsun(p) = 0.;
vegetation.mlcanopyinst.etvegsha(p) = 0.;
vegetation.mlcanopyinst.gppvegsun(p) = 0.;
vegetation.mlcanopyinst.gppvegsha(p) = 0.;

for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
    if (vegetation.canopy.dpai(p,ic) > 0)
        vegetation.mlcanopyinst.irvegsun(p) = vegetation.mlcanopyinst.irvegsun(p) + vegetation.mlcanopyinst.irleaf(p,ic) * vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.irvegsha(p) = vegetation.mlcanopyinst.irvegsha(p) + vegetation.mlcanopyinst.irleaf(p,ic) * vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.shvegsun(p) = vegetation.mlcanopyinst.shvegsun(p) + vegetation.mlcanopyinst.shleaf(p,ic,isun) * vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.shvegsha(p) = vegetation.mlcanopyinst.shvegsha(p) + vegetation.mlcanopyinst.shleaf(p,ic,isha) * vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.lhvegsun(p) = vegetation.mlcanopyinst.lhvegsun(p) + vegetation.mlcanopyinst.lhleaf(p,ic,isun) * vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.lhvegsha(p) = vegetation.mlcanopyinst.lhvegsha(p) + vegetation.mlcanopyinst.lhleaf(p,ic,isha) * vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.etvegsun(p) = vegetation.mlcanopyinst.etvegsun(p) + (vegetation.mlcanopyinst.evleaf(p,ic,isun) + vegetation.mlcanopyinst.trleaf(p,ic,isun)) * vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic);
        vegetation.mlcanopyinst.etvegsha(p) = vegetation.mlcanopyinst.etvegsha(p) + (vegetation.mlcanopyinst.evleaf(p,ic,isha) + vegetation.mlcanopyinst.trleaf(p,ic,isha)) * vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic);
        fracgreen = vegetation.mlcanopyinst.fdry(p,ic) / (1  - vegetation.mlcanopyinst.fwet(p,ic));
        vegetation.mlcanopyinst.gppvegsun(p) = vegetation.mlcanopyinst.gppvegsun(p) + vegetation.mlcanopyinst.ag(p,ic,isun) * vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic) * fracgreen;
        vegetation.mlcanopyinst.gppvegsha(p) = vegetation.mlcanopyinst.gppvegsha(p) + vegetation.mlcanopyinst.ag(p,ic,isha) * vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic) * fracgreen;
    end
end
end

