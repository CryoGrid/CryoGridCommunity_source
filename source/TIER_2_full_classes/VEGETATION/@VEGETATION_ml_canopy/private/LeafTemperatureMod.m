function [vegetation] = LeafTemperatureMod (vegetation, p, ic, il)
%p, ic, il, shr_kind, abortutils, clm_varctl, WaterVaporMod)
%-----------------------------------------------------------------------
% %DESCRIPTION:
% Leaf temperature and energy fluxes (excluding evaporation of intercepted water)

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Leaf Temperature DONE
%Leaf Temperature and Energy Fluxes
%Input
% vegetation.mlcanopyinst.tleaf = vegetation.mlcanopyinst.tleaf;
% dpai = vegetation.canopy.dpai;                                  % vegetation.mlcanopyinst%dpai         % Layer plant area index (m2/m2)
% vegetation.mlcanopyinst.cpair = vegetation.mlcanopyinst.cpair;                                % vegetation.mlcanopyinst%cpair        % Specific heat of air at constant pressure, at reference height (J/mol/K)

%eref instead of vegetation.mlcanopyinst.eair
% eair = vegetation.mlcanopyinst.eair;

% gs = vegetation.mlcanopyinst.gs;
% rnleaf = vegetation.mlcanopyinst.rnleaf;
% tleaf_old = vegetation.mlcanopyinst.tleaf_old;%_old;
% cpleaf = vegetation.mlcanopyinst.cpleaf;
% fwet = vegetation.mlcanopyinst.fwet ;        % Fraction of plant area index that is wet
% fdry = vegetation.mlcanopyinst.fdry ;        % Fraction of plant area index that is green and dry
% gbh  = vegetation.mlcanopyinst.gbh  ;        % Leaf boundary layer conductance, heat (mol/m2 leaf/s)
% gbv  = vegetation.mlcanopyinst.gbv ;         % Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)gs = 1;

%Output                                                         % die kommen als output raus. so au zeile 696 in scaler profile schreiben.
% stleaf = vegetation.mlcanopyinst.stleaf;
% shleaf = vegetation.mlcanopyinst.shleaf;
% lhleaf = vegetation.mlcanopyinst.lhleaf;

% evleaf = vegetation.mlcanopyinst.evleaf;       % Leaf evaporation flux (mol H2O/m2 leaf/s)
% trleaf = vegetation.mlcanopyinst.trleaf;       % Leaf transpiration flux (mol H2O/m2 leaf/s)
% tref  = vegetation.mlcanopyinst.tref;         % Air temperature at reference height (K)
% pref  = vegetation.mlcanopyinst.pref;         % Air pressure at reference height (Pa)
% tair  = vegetation.mlcanopyinst.tair;         % Air temperature profile (K)
% tleaf  = vegetation.mlcanopyinst.tleaf;        % Leaf temperature (K)

% Time step (s)
dtime = vegetation.params.dtime_sub;

if (vegetation.canopy.dpai(p,ic) > 0) % >0  leaf layer
    
    %Saturation vapor pressure (Pa -> mol/mol)
    [esat, desat] = Satvap (vegetation.mlcanopyinst.tleaf_old(p,ic,il));
    
    qsat = esat ./ vegetation.mlcanopyinst.pref(p);
    dqsat = desat ./ vegetation.mlcanopyinst.pref(p);
    
    % Latent heat of vaporization
    
    [lambda] = LatVap(vegetation.mlcanopyinst.tref(p), vegetation);
    
    
    % Leaf conductance for transpiration
    
    gleaf = vegetation.mlcanopyinst.gs(p,ic,il) .* vegetation.mlcanopyinst.gbv(p,ic,il) ./ (vegetation.mlcanopyinst.gs(p,ic,il) + vegetation.mlcanopyinst.gbv(p,ic,il));
    
    % Total conductance (including evaporation)
    
    gw = gleaf .* vegetation.mlcanopyinst.fdry(p,ic) + vegetation.mlcanopyinst.gbv(p,ic,il) .* vegetation.mlcanopyinst.fwet(p,ic);
   
    % Linearized leaf temperature calculation that balances the energy budget    
    num1 = 2 .* vegetation.mlcanopyinst.cpair(p) .* vegetation.mlcanopyinst.gbh(p,ic,il); %1.3; % in book vegetation.mlcanopyinst.gbh(p,ic,il) = 1.3, here it was only 0.0012
    num2 = lambda .* gw;
    num3 = vegetation.mlcanopyinst.rnleaf(p,ic,il) - lambda .* gw .* (qsat - dqsat .* vegetation.mlcanopyinst.tleaf_old(p,ic,il)) ...
        + vegetation.mlcanopyinst.cpleaf(p,ic) ./ dtime .* vegetation.mlcanopyinst.tleaf_old(p,ic,il);      %vegetation.mlcanopyinst.cpleaf = heat capacity
    den = vegetation.mlcanopyinst.cpleaf(p,ic) ./ dtime + num1 + num2 .* dqsat;
    vegetation.mlcanopyinst.tleaf(p,ic,il) = (num1 .* vegetation.mlcanopyinst.tair(p,ic) + num2 .* vegetation.mlcanopyinst.eair(p,ic) ./ vegetation.mlcanopyinst.pref(p) + num3) ./ den;
    
    if vegetation.mlcanopyinst.tleaf(p,ic,il)<0
     disp('shit')
    end
    
    % Storage flux
    vegetation.mlcanopyinst.stleaf(p,ic,il) = (vegetation.mlcanopyinst.tleaf(p,ic,il) - vegetation.mlcanopyinst.tleaf_old(p,ic,il)) .* vegetation.mlcanopyinst.cpleaf(p,ic) ./ dtime;

    % Sensible heat flux
    vegetation.mlcanopyinst.shleaf(p,ic,il) = 2 .* vegetation.mlcanopyinst.cpair(p) .* (vegetation.mlcanopyinst.tleaf(p,ic,il) - vegetation.mlcanopyinst.tair(p,ic)) .* vegetation.mlcanopyinst.gbh(p,ic,il);
    
    % Transpiration and evaporation water fluxes: mol H2O/m2/s
    
    num1 = qsat + dqsat .* (vegetation.mlcanopyinst.tleaf(p,ic,il) - vegetation.mlcanopyinst.tleaf_old(p,ic,il)) - vegetation.mlcanopyinst.eair(p,ic) ./ vegetation.mlcanopyinst.pref(p);
    vegetation.mlcanopyinst.trleaf(p,ic,il) = gleaf .* vegetation.mlcanopyinst.fdry(p,ic) .* num1;
    vegetation.mlcanopyinst.evleaf(p,ic,il) = vegetation.mlcanopyinst.gbv(p,ic,il) .* vegetation.mlcanopyinst.fwet(p,ic) .* num1;

    
    % Latent heat flux
    % ev = leaf evaporation flux
    % tr = transpiration flux
    
    vegetation.mlcanopyinst.lhleaf(p,ic,il) = (vegetation.mlcanopyinst.trleaf(p,ic,il) + vegetation.mlcanopyinst.evleaf(p,ic,il)) .* lambda;
       
    % Error check
    
    % leaf net radiation - sensible heat flux - latent heat flux - storage flux
%     if max(vegetation.mlcanopyinst.lhleaf)>20
%         a=1
%     end
    err = vegetation.mlcanopyinst.rnleaf(p,ic,il) - vegetation.mlcanopyinst.shleaf(p,ic,il) - vegetation.mlcanopyinst.lhleaf(p,ic,il) - vegetation.mlcanopyinst.stleaf(p,ic,il);
    if (abs(err) > 1e-3)
        %         disp(ic);
        disp('LeafTemperatureMod error:');
        disp(err);
        %         write (iulog,*) vegetation.mlcanopyinst.rnleaf(p,ic,il), shleaf(p,ic,il), lhleaf(p,ic,il), stleaf(p,ic,il), err
        %         write (iulog,*) tleaf(p,ic,il), vegetation.mlcanopyinst.tleaf_old(p,ic,il)
        %         disp('ERROR: LeafTemperatureMod: energy balance error');
    end
    
else % non-leaf layer
    
    vegetation.mlcanopyinst.tleaf(p,ic,il) = vegetation.mlcanopyinst.tair(p,ic);
    %CHANGE SEBAS
    %vegetation.mlcanopyinst.tleaf(p,1,il) = NaN; %vegetation.mlcanopyinst.tg(p);  %SS: tleaf(p,1) würde eigentlich gar nicht existieren.. 
    vegetation.mlcanopyinst.tleaf(p,1,il) = vegetation.mlcanopyinst.tair(p,1);
    %end change SEBAS
    vegetation.mlcanopyinst.stleaf(p,ic,il) = 0.;
    vegetation.mlcanopyinst.shleaf(p,ic,il) = 0.;
    vegetation.mlcanopyinst.lhleaf(p,ic,il) = 0.;
    vegetation.mlcanopyinst.evleaf(p,ic,il) = 0.;
    vegetation.mlcanopyinst.trleaf(p,ic,il) = 0.;
    
end
end