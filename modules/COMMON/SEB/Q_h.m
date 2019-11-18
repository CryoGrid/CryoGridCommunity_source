function Q_h = Q_h(ground, forcing)
%q=specific humidity [kg/kg]

uz = forcing.TEMP.wind;
z =  ground.PARA.airT_height; 
z0 = ground.PARA.z0;
Tz = forcing.TEMP.Tair;
TForcing = ground.STATVAR.T(1);
Lstar = ground.STATVAR.Lstar; 
p = forcing.TEMP.p;


Tz=Tz+273.15;
TForcing=TForcing+273.15;

rho = p./(287.058.*Tz); %air density [kg m^(-3)]
cp=1005;
kappa=0.4;
g=9.81;
sigma=5.67e-8;


Q_h  =-rho.*cp.*kappa.* uz.*kappa./(log(z./z0)- psi_M(ground, z./Lstar, z0./Lstar)) .* (Tz-TForcing)./(log(z./z0)- psi_H(ground, z./Lstar, z0./Lstar));