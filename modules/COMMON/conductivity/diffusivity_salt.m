function ground = diffusivity_salt(ground)

water = ground.STATVAR.water./ground.STATVAR.layerThick;

D0 = ((6.06 + 9.60)/2  + max(ground.STATVAR.T, 0) .* (0.297  + 0.438)/2) .* 1e-10; %from Boudreau, B., 1997, Diagenetic Models and thier implementation, Springer, Berlin.
%average between values for Na+ and Cl-
tortuosity = ground.PARA.tortuosity; % some nice value

ground.STATVAR.diffusivitySalt = D0 .* water ./tortuosity.^2;