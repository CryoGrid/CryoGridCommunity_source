function [vegetation] = SoilFluxesMultilayer (vegetation, p)

c = p;
ivis = vegetation.params.vis; % Array index for visible waveband
inir = vegetation.params.nir; % Array index for near-infrared waveband

%  Current ground temperature

% Simone
% if vegetation.mlcanopyinst.snow == 1
% vegetation.mlcanopyinst.tg(p) = vegetation.mlcanopyinst.tg_snow(1); %p,0
% else
%SEBAS:removed
%vegetation.mlcanopyinst.tg(p) = vegetation.mlcanopyinst.tair_old(p,1); %p,0
%end SEBAS:removed
% end

%  Net radiation

vegetation.mlcanopyinst.rnsoi(p) = vegetation.flux.swsoi(p,ivis) + vegetation.flux.swsoi(p,inir) + vegetation.flux.irsoi(p);

%  Latent heat of vaporization

[lambda] = LatVap(vegetation.mlcanopyinst.tref(p), vegetation);

%  Relative humidity in soil airspace
%  Relative humidity of airspace at soil surface (fraction) = exp(Gravitational acceleration (m/s2) * Molecular mass of water (kg/mol) * Soil layer matric potential (mm) / (universal gas constant
%  [J/K/kmole] * Soil temperature (K))
%  fraction = exp(m/s2) * kg/mol * mm / (J/K/kmole * K) 
%  vegetation.mlcanopyinst.rhg(p) = exp(vegetation.physcon.grav *  vegetation.physcon.mmh2o *  smp_l(c,1)*1.e-03 / ( vegetation.physcon.rgasc * vegetation.physcon.t_soisno(c,1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserdampfdruck über schnee muss noch geändert werden.
% molar_volume = vegetation.physcon.rgasc * vegetation.soilvar.t_soisno(c,1) / 101.3; % [m3/mol] 
%molar_volume = vegetation.physcon.rgasc * vegetation.soilvar.t_top_surfacecell / 101.3; % [m3/mol] 
% vegetation.mlcanopyinst.rhg(p) = exp(vegetation.soilvar.soil_water_matric_potential(c,1) * molar_volume / (vegetation.physcon.rgasc * vegetation.soilvar.t_soisno(c,1)));
%vegetation.mlcanopyinst.rhg(p) = exp(vegetation.soilvar.soil_water_matric_potential(c,1) * molar_volume / (vegetation.physcon.rgasc * vegetation.soilvar.t_top_surfacecell)); %vegetation.soilvar.t_soisno(c,1)));

%SEBAS changed
%The lines above are based on CLM 5 formulations, CLM 4.5 is different,
%compare Eq. 5.73 in CLM 4.5 manual - I have chenged to CLM 4.5
vegetation.mlcanopyinst.rhg(p) = exp(vegetation.soilvar.soil_water_matric_potential(c,1) .* 9.81 ./ ((8.3145e+03 ./ 18.016) .*vegetation.soilvar.t_top_surfacecell));
%added beta factor from CLM 4.5, this is really what seems to supress E
beta = 1 +  double(vegetation.soilvar.h2osoi_vol(c,1)<0.5) .* (-1 +  0.25 .* (1-(cos(pi() .* vegetation.soilvar.h2osoi_vol(c,1) ./ 0.5))).^2);
%end SEB changed

% % rgasc back to old value of  8.3145e+03: --> doesn´t matter, molar volume * e+03 / rgasce+03 = same as without e+03
% molar_volume = 8.3145e+03 * vegetation.soilvar.t_soisno(c,1) / 101.3; % [m3/mol] 
% vegetation.mlcanopyinst.rhg(p) = exp(vegetation.soilvar.soil_water_matric_potential * molar_volume / (8.3145e+03 * vegetation.soilvar.t_soisno(c,1)));

%  Soil conductance to water vapour diffusion
gws = 1. / vegetation.mlcanopyinst.soilresis(c);                      % ! s/m -> m/s                        %this is a constant - why?
gws = gws * vegetation.mlcanopyinst.rhomol(p);                         %   ! m/s -> mol H2O/m2/s            %does not change too much either
gw = vegetation.mlcanopyinst.ga_prof(p,1) * gws / (vegetation.mlcanopyinst.ga_prof(p,1) + gws);  % ! total conductance
%gw is called g0 in scalar_fluxes!! But it seems to be exactly the same
%THIS IST STILL WRONG!! Contains some of CLM5 formulation, now not used

%  Saturation vapor pressure at ground temperature (Pa -> mol/mol)
[esat, desat] = Satvap(vegetation.mlcanopyinst.tg(p));
qsat = esat / vegetation.mlcanopyinst.pref(p) ;
dqsat = desat / vegetation.mlcanopyinst.pref(p);

%SEBAS:removed
%  Calculate soil surface temperature
% num1 = vegetation.mlcanopyinst.cpair(p) * vegetation.mlcanopyinst.ga_prof(p,1); %p,0
% num2 = lambda * gw;
% % num3 = vegetation.soilvar.thk(c) / (vegetation.soilvar.z(c)-vegetation.soilvar.zi(c)); %no snow layer, so leave out snl(c)+1
% num3 = vegetation.soilvar.thk_topsurfacecell / (vegetation.soilvar.dz_topsurfacecell); %no snow layer, so leave out snl(c)+1
% % num4 = vegetation.mlcanopyinst.rnsoi(p) - num2 * vegetation.mlcanopyinst.rhg(p) * (qsat - dqsat * vegetation.mlcanopyinst.tg(p)) + num3 * vegetation.soilvar.t_soisno(c,1);
% num4 = vegetation.mlcanopyinst.rnsoi(p) - num2 * vegetation.mlcanopyinst.rhg(p) * (qsat - dqsat * vegetation.mlcanopyinst.tg(p)) + num3 * vegetation.soilvar.t_top_surfacecell;
% 
% den = num1 + num2 * dqsat * vegetation.mlcanopyinst.rhg(p) + num3;
% vegetation.mlcanopyinst.tg(p) = (num1*vegetation.mlcanopyinst.tair(p,2) + num2*vegetation.mlcanopyinst.eair(p,2)/vegetation.mlcanopyinst.pref(p) + num4) / den; %p,1
%end SEBAS:removed


%  Sensible heat flux
vegetation.mlcanopyinst.shsoi(p) = vegetation.mlcanopyinst.cpair(p) * (vegetation.mlcanopyinst.tg(p) - vegetation.mlcanopyinst.tair(p,2)) * vegetation.mlcanopyinst.ga_prof(p,1); %p,1%p,0

%SEBAS:removed
%  Latent heat flux - remember that tair_old(p,0) is tg(p) before the update for time n+1
%vegetation.mlcanopyinst.eg(p) = vegetation.mlcanopyinst.rhg(p) * (esat + desat * (vegetation.mlcanopyinst.tg(p) - vegetation.mlcanopyinst.tair_old(p,1))); %p,0
%end SEBAS:removed

vegetation.mlcanopyinst.eg(p) = vegetation.mlcanopyinst.rhg(p) * esat;
vegetation.mlcanopyinst.lhsoi(p) = lambda / vegetation.mlcanopyinst.pref(p) .* beta .*(vegetation.mlcanopyinst.eg(p) - vegetation.mlcanopyinst.eair(p,2)) * gw; %p,1

%  Soil heat flux
% vegetation.mlcanopyinst.gsoi(p) = vegetation.soilvar.thk(c) * (vegetation.mlcanopyinst.tg(p) - vegetation.soilvar.t_soisno(c)) / (vegetation.soilvar.z(c)-vegetation.soilvar.zi(c));
% %SEBAS:removed
% vegetation.mlcanopyinst.gsoi(p) = vegetation.soilvar.thk_topsurfacecell(c) * (vegetation.mlcanopyinst.tg(p) - vegetation.soilvar.t_top_surfacecell) / (vegetation.soilvar.dz_topsurfacecell);
% 
% %  Error check
% 
% err = vegetation.mlcanopyinst.rnsoi(p) - vegetation.mlcanopyinst.shsoi(p) - vegetation.mlcanopyinst.lhsoi(p) - vegetation.mlcanopyinst.gsoi(p);
% if (abs(err) > 0.001)
% %     disp(' ERROR: SoilFluxesMultilayerMod: energy balance error');
% end
%end SEBAS:removed

%  Water vapor flux: W/m2 -> mol H2O/m2/s
vegetation.mlcanopyinst.etsoi(p) = vegetation.mlcanopyinst.lhsoi(p) / lambda;

end