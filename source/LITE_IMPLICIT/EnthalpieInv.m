function [T, fl] = EnthalpieInv(H, WaterIce, HC)

rho = 1000; %[kg/m^3]
Lsl = 334000; %[J/kg]
L = rho*Lsl;%[J/m^3]

theta = WaterIce;
%avoid zeros division
theta(theta==0)=1e-8;

fl = (H>0 & H<=L*theta).*(H./(L*theta)) + (H>L*theta);
%c = theta.*(fl.*cl + (1-fl).*cs) + (1-theta).*cm;
c = HC;

T = (H>=(L*theta)).*((H-L*theta)./(c)) + (H<0).*(H./(c));

end
