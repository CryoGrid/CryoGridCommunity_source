function [H, dHdT] = Enthalpie(T, Water, HC)

rho = 1000; %[kg/m^3]
Lsl = 334000; %[J/kg]
L = rho*Lsl;%[J/m^3]

theta = Water; %[Vol. fraction]

H = (T.*HC + theta.*L); %[J/m^3]

%dHdT = zeros(size(T));
dHdT = HC.*(T~=0.0) + 1e8.*(T==0.0); %[J/m^3K]
%dHdT = HC; %[J/m^3K]
end
