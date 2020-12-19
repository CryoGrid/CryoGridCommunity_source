
theta_air =[0.01:0.01:0.4]';

porosity = 0.85;
T = 276;
R = 8.31446261815324;
g = 9.81;
d_theta = 0.001;
p0 = 1000e2 + 1.5.*1000.*g; 

p0/1000/g

theta_water = porosity - theta_air;


res1 = (waterPotential(theta_water+d_theta, porosity) - waterPotential(theta_water - d_theta, porosity) ) ./ (2.*d_theta) .* R .* T .* theta_air ./ (waterPotential(theta_water, porosity)+p0/1000/g);
res2 = R.*T;

plot(theta_air, [res1+res2])

function waterPot = waterPotential(waterIce, porosity)

alpha = 2.31;
n = 1.29;

m=1-1./n;
waterPot = -1./alpha .*(((waterIce)./(porosity)).^(-1./m) -1).^(1./n);

end