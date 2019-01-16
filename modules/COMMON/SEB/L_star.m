function ground = L_star(ground, forcing)


uz = forcing.TEMP.wind;
z =  ground.PARA.airT_height; 
z0 = ground.PARA.z0;
Tz = forcing.TEMP.Tair+273.15;
Lstar = ground.STATVAR.Lstar; 
p = forcing.TEMP.p;
Qh = ground.STATVAR.Qh;
Qe = ground.STATVAR.Qe;

  

    
%rho=1.293;
rho = p./(287.058.*Tz); %air density [kg m^(-3)]
cp=1005;
%L=2.8*10^6;  %changed
if Tz >=273.15
    L=1e3.*(2500.8 - 2.36.*(Tz-273.15));  %latent heat of evaporation of water
else
    L=1e3.*2834.1; %latent heat of sublimation
end


kappa=0.4;
g=9.81;


u_star = real(uz.*kappa./(log(z./z0)- psi_M(ground, z./Lstar, z0./Lstar)));
L_star = real(-rho.*cp.*Tz./kappa./g.*u_star.^3./(Qh + 0.61.*cp./L.*Tz.*Qe));
L_star=(abs(L_star)<1e-7).*L_star./abs(L_star).*1e-7 + (abs(L_star)>=1e-7).*L_star;  %limits Lstar

ground.STATVAR.Lstar = L_star;
ground.STATVAR.u_star = u_star;



% if isnan(L_star)
%     jklkl
% end



