function visc = viscosity_water(T)

T = max(0,T);
T= T+273.15;

%from Wikepedia
A = 1.856e-11 .* 1e-3;
B = 4209;
C = 0.04527;
D = -3.376e-5;


visc = A.*exp(B./T + C.* T + D.* T.^2); %[Pa sec]
