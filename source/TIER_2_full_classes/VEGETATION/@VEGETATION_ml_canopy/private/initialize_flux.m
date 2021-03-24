function [ground] = initialize_flux(ground)

% nleaf = ground.STATVAR.vegetation.mlcanopyinst.nleaf;
% ncan = ground.STATVAR.vegetation.mlcanopyinst.ncan;
ground.STATVAR.vegetation.mlcanopyinst.numrad = 2;                             % Number of wavebands
numrad = ground.STATVAR.vegetation.mlcanopyinst.numrad;

% ground.STATVAR.vegetation.flux.albsoib = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
% ground.STATVAR.vegetation.flux.albsoid = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.flux.fracsun = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan);
ground.STATVAR.vegetation.flux.fracsha = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan);
ground.STATVAR.vegetation.flux.swsoi = zeros(1,2);
ground.STATVAR.vegetation.flux.swveg = zeros(1,2);
ground.STATVAR.vegetation.flux.swvegsun = zeros(1,2);
ground.STATVAR.vegetation.flux.swvegsha = zeros(1,2);
ground.STATVAR.vegetation.flux.swleaf = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf,numrad);
ground.STATVAR.vegetation.flux.albcan = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.flux.apar = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.flux.ir_source_sun = zeros(1);
ground.STATVAR.vegetation.flux.ir_source_sha = zeros(1);
ground.STATVAR.vegetation.flux.ir_source = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.flux.irsoi = zeros(1);
ground.STATVAR.vegetation.flux.irveg = zeros(1);
ground.STATVAR.vegetation.flux.irup = zeros(1);

end 