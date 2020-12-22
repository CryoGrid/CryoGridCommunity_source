theta_m= 0;

for theta_o = [0.05:0.05:0.95]
 
    porosity = 1 - theta_o - theta_m;

    c_w = 4.2e6;
    c_i = 1.9e6;
    c_m = 2e6;
    c_o = 2.5e6;
    L_sl = 3.34e8;
    T0=273.15;
    beta_interface = 1./3;

alpha = 8e-4; %(Pa^-1)
m=0.19;


 alpha = 1.11e-4;
 n=1.48;
 alpha = 1.49e-4;
 n=1.25;
 m=1-1./n;

% sat_waterIce_list = [0.2; 0.5; 1];
 T=[-40:0.01:0]';


n=1./(1-m);

% for i=1:3
%     sat_waterIce = sat_waterIce_list(i)
sat_waterIce_list = [0.01:0.01:1];
 for i=1:size(sat_waterIce_list,2)
    
sat_waterIce = sat_waterIce_list(i);
% sat_waterIce=1;

mwp0 = 1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n);

A = -L_sl.*T./T0 .*beta_interface .* double(T<0) + mwp0;

sat_water = double(A>0) .* (1+(alpha.*A).^n).^(-m) + double(A<=0);
sat_ice = sat_waterIce-sat_water;
% plot(T,sat_water)
% hold on

energy = T.* ((1-porosity) .* 2e6 + porosity .* sat_waterIce .* (c_w .* double(T>=0)+ c_i.* double(T<0)));
energy = energy - double(T<0) .* L_sl .* porosity .* (sat_waterIce-sat_water);

%   plot(energy, T)
  
  
C0 = (1-porosity).* 2e6 + porosity .* sat_waterIce .* c_i;
X = -L_sl./ T0 .* beta_interface;
L0 = porosity.* L_sl;

E_prime = energy ./ L0 + sat_waterIce + C0 .* mwp0 ./ X./ L0;
T_prime = alpha .* (X .* T + mwp0);

gamma = C0 ./ alpha ./ X ./ L0 ;

T_prime2 = (0:2000)';
E_prime2 = gamma .* T_prime2 + (1 + T_prime2.^n).^(-m);
% plot(E_prime2, T_prime2, '+')
% 
% % figure
%   hold on
%   plot(E_prime, T_prime)
plot(sat_waterIce, log(-gamma), '+')
hold on
 end
end

porosity_max = 0.95;
gamma_max =  ((1-porosity_max).* c_o ) ./ alpha ./ X ./ (porosity_max.* L_sl);
log(-gamma_max)
porosity_min = 0.05;
gamma_min =  ((1-porosity_min).* c_m  + c_i) ./ alpha ./ X ./ (porosity_min.* L_sl);
log(-gamma_min)



% end
  
%  plot(energy, sat_water)



% B = sat_water./(1-sat_ice);
% P_cgl = 1./alpha .* ((B.^(-1./m)-1)).^(1./n);
% P_csl = -L_sl.*T./T0 .*double(T<0);

% get_E(-10)
% 
% 
% 
% fsolve(@get_E, -15)

function energy = get_E(T)

porosity = 0.5;
alpha = 8e-4; %(Pa^-1)
m=0.19;
L_sl = 3.34e8;
T0=273.15;

alpha = 1.11e-4;
n=1.48;
m=1-1./n;

% sat_waterIce_list = [0.2; 0.5; 1];


n=1./(1-m);

% for i=1:3
%     sat_waterIce = sat_waterIce_list(i)

sat_waterIce = 0.5;
A = -L_sl.*T./T0 ./3 .*double(T<0) + 1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n);

sat_water = double(A>0) .* (1+(alpha.*A).^n).^(-m) + double(A<=0);
sat_ice = sat_waterIce-sat_water;


energy = T.* (porosity .* 2e6 + (1-porosity) .* sat_waterIce .* (4.2e6 .* double(T>=0)+ 1.9e6.* double(T<0)));
energy = energy - double(T<0) .* 3.34e8 .* (1-porosity) .* sat_ice;

energy= energy +8.9414e+07;
end