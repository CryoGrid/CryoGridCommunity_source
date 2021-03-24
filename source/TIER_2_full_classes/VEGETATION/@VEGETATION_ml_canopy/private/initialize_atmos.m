function [ground] = initialize_atmos(ground)

ground.STATVAR.vegetation.mlcanopyinst.numrad = 2;                             % Number of wavebands
numrad = ground.STATVAR.vegetation.mlcanopyinst.numrad;

ground.STATVAR.vegetation.atmos.swskyb = zeros(1,numrad);
ground.STATVAR.vegetation.atmos.swskyd = zeros(1,numrad);
end