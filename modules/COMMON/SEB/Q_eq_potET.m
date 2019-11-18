function Q_e = Q_eq_potET(ground, forcing)


uz = forcing.TEMP.wind;
p = forcing.TEMP.p;
q = forcing.TEMP.q;
Tz = forcing.TEMP.Tair;


z =  ground.PARA.airT_height; 
z0 = ground.PARA.z0;

TForcing = ground.STATVAR.T(1);
Lstar = ground.STATVAR.Lstar; 





%q = specific humidity [kg/kg]
Tz=Tz+273.15;
TForcing=TForcing+273.15;
%p=1005*100; %air preassure [Pa] 
%rho=1.293;
%rho = (p-(p.*q))./(287.058.*Tz) + (p.*q)./(461.495.*Tz); %air density [kg m^(-3)]
rho = p./(287.058.*Tz); %air density [kg m^(-3)]
cp=1005;
L_w=1e3.*(2500.8 - 2.36.*(TForcing-273.15));  %latent heat of evaporation of water
L_i=1e3.*2834.1; %latent heat of sublimation
kappa=0.4;
g=9.81;




if TForcing<273.15
    Q_e = -rho.*L_i.*kappa.*uz.*kappa./(log(z./z0)- psi_M(ground, z./Lstar, z0./Lstar)).*(q-satPresIce(ground, TForcing)./p)./(log(z./z0)- psi_H(ground, z./Lstar, z0./Lstar));
else
    Q_e = -rho.*L_w.*kappa.*uz.*kappa./(log(z./z0)- psi_M(ground,z./Lstar, z0./Lstar)).*(q-satPresWater(ground, TForcing)./p)./(log(z./z0)- psi_H(ground, z./Lstar, z0./Lstar));
end
