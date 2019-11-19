function ground = finalize_STATVAR(ground)


ground.STATVAR.Lstar = -100;
ground.STATVAR.Qh = 0;
ground.STATVAR.Qe = 0;


%%NC added - SEB initializers

size_L_star_smoothing=1;

ground.STATVAR.Qnet=0;
ground.STATVAR.Qg=0;
ground.STATVAR.Sout=0;
ground.STATVAR.Lout=0;
ground.PARA.newSnow=0;
ground.PARA.SEB.L_star=-100000+zeros(1,size_L_star_smoothing);
ground.PARA.u_star=10;


ground.PARA.Qsurf = 0;  % for EB checks
