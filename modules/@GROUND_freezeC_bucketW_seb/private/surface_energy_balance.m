function ground = surface_energy_balance(ground, forcing)


ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
ground.STATVAR.Qh = Q_h(ground, forcing);
ground.STATVAR.Qe = Q_eq(ground, forcing);

ground.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;