
function [vegetation] = CanopyFluxesSum (vegetation, p) %, ic, il
%
% %DESCRIPTION:
% Sum leaf and soil fluxes

ivis = vegetation.params.vis; % Array index for visible waveband
inir = vegetation.params.nir; % Array index for near-infrared waveband
isun = vegetation.params.sun; % Array index for sunlit leaf
isha = vegetation.params.sha; % Array index for shaded leaf


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

avail = radin - radout - vegetation.mlcanopyinst.gsoi(p);
flux = vegetation.mlcanopyinst.shflx(p) + vegetation.mlcanopyinst.lhflx(p) + vegetation.mlcanopyinst.stflx(p);
err = avail - flux;
if (abs(err) > 0.01)
    disp ('ERROR: CanopyFluxesSum: energy conservation error (3)');
    disp (err);
end

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

