function  ground = surfaceEnergyBalanceInfiltration(ground, forcing)

% T = ground.STATVAR.T;
% wc = ground.STATVAR.water ./ ground.STATVAR.layerThick;
% Lstar=ground.STATVAR.Lstar;
% z = forcing.PARA.airT_height;

L = 1e3.*(2500.8 - 2.36.*ground.STATVAR.T(1)) .* ground.CONST.rho_w;  %latent heat of evaporation of water



ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
ground.STATVAR.Qh = Q_h(ground, forcing);
Q_potET = Q_eq_potET(ground, forcing);

% 2 cases  for weighting factor
%-> 1. first cell frozen OR Q_eq_potET <0 (condensation) - > use Q_eq_potET
%-> else do the depth-dependent weighting factor stuff

if ground.STATVAR.T(1) <= 0 || Q_potET <= 0
    weighting_factor_ET = 1;
else
    fraction_T = getT_fraction(ground);
    fraction_E = getE_fraction(ground);
    
    depth_weighting_E = exp(-1./ground.PARA.evaporationDepth.* (cumsum(ground.STATVAR.layerThick) - ground.STATVAR.layerThick(1)./2));  %exponential decrease with depth ,needs to be corrected for subsieded grid
    depth_weighting_E(depth_weighting_E<0.05)=0;

    depth_weighting_E=depth_weighting_E.*ground.STATVAR.layerThick./sum(depth_weighting_E.*ground.STATVAR.layerThick); %normalize
    
    depth_weighting_T=exp(-1./ground.PARA.rootDepth.*(cumsum(ground.STATVAR.layerThick) - ground.STATVAR.layerThick(1)./2)); 
    depth_weighting_T(depth_weighting_T<0.05)=0;

    depth_weighting_T=depth_weighting_T .* ground.STATVAR.layerThick./sum(depth_weighting_T.*ground.STATVAR.layerThick);
    
    fraction_ET = fraction_T.*depth_weighting_T.*ground.PARA.ratioET + fraction_E.*depth_weighting_E.*(1-ground.PARA.ratioET);
    
    weighting_factor_ET = sum(fraction_ET, 1);
    
end

ground.STATVAR.Qe = weighting_factor_ET .* Q_potET;

%energy balance
ground.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;



%then do the water calculation

ground.TEMP.dwc_dt = ground.STATVAR.T .* 0;
if ground.STATVAR.T(1) >= 0
    if ground.STATVAR.Qe > 0
    
        fraction_ET=fraction_ET./sum(fraction_ET, 1); %normalize
        
        ground.TEMP.dwc_dt = -ground.STATVAR.Qe./L.*fraction_ET;    %in m water per sec
    else  %condensation
        ground.TEMP.dwc_dt(1) = -ground.STATVAR.Qe ./ L;
    end
end
    
    
    
