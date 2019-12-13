function [liqWater, Tmelt] = watercontent_temperature_salt(ground, optionalArgs)
%T is required to be in Kelvin!
%convert to Kelvin, if you use ground.STATVAR.T or 
%give T in Kelvin via optionalArgs

if nargin > 1 %optionalArgs
    %T is in Kelvin
    T = optionalArgs{1};
    saltConc = optionalArgs{2};
else
    %Convert from degreeC to Kelvin
    T = ground.STATVAR.T + 273.15;      % temperature
    saltConc = ground.STATVAR.saltConc; % salt concentration
end

theta_sat = ground.STATVAR.water;   % total water content

%see del'amico
alpha = ground.PARA.alpha;         %in Van genuchten model [m^-1]
n = ground.PARA.n;                           %in van genuchten model [-]


%we already did this, I think? ns = 2*salinity?
%ns=2.*ground.STATVAR.ns; %one mole NaCl forms 2 moles of ions

R = ground.CONST.R;
g = ground.CONST.g;
rho_w = ground.CONST.rho_w;
L_f = ground.CONST.L_f;
Tmelt_free_water = ground.CONST.Tmelt_free_water; % in K

Tmelt = Tmelt_free_water + Tmelt_free_water./L_f.*(-R.*saltConc.*Tmelt_free_water); %in Kelvin

%Eq. 17 from Dall Amico, but using the solute potential of the UNFROZEN
%soil as the last term - I think this is correct now consistent with
%freezing_point_depression in free water

%L_f./(rho_w.*g.*Tmelt).*(T-Tmelt).*(T < Tmelt);
%R.*ns.*T./(rho_w.*g.*liqWater)-R.*ns.*T./(rho_w.*g.*theta_sat);

water_pot = L_f./(rho_w.*g.*Tmelt).*(T-Tmelt).*(T < Tmelt);  %+R.*ns.*(T-Tmelt)./(rho_w.*g).*(T < Tmelt);
% the first term is equivalentto Eq. 2 in Dall Amico - this is the effect of the matric potential
% the third term is the solute potential of the unfrozen soil - somehow this is needed as a reference potential
% the sethermCond term is the solute potential of the freezing soil - this is where the magig happens 
% as it increases with decreasing liqWater - but not sure why this is PLUS -
% but it must be plus to offset the efect of the first term which is minus
% and becoming more negative with decreasing T - otherwise the freezing
% happens even faster

water_pot = water_pot.*(water_pot<0);
% this is needed to get unchanged water contents above Tmelt

liqWater = theta_sat./(1+(alpha.*abs(water_pot)).^n).^(1-1./n);
%this is the final formula, but waterPot depends also on liqWater, so it's a
%coupled eq. system

%res=liqWater-liqWater_modeled;
% this is needed for fsolve, if res=0, the equation system is solved

end

