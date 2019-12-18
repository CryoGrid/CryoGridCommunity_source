function ground = add_water_GROUND(ground,waterflux)
% SAME AS BucketScheme.m used by SNOW_freezeC_(..), based on CryoGrid3 (J. Nitzbonn ?)
% R. B. Zweigel, December 2019

T = ground.STATVAR.T; 
wc = ground.STATVAR.water ./ ground.STATVAR.layerThick;
dwc = ground.STATVAR.water.*0;  %should be in m
dwc(1) = waterflux;

K_delta = ground.STATVAR.layerThick;  %in m
porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick; 

fieldCapacity = ground.STATVAR.field_capacity; 



i=1;

while  T(i)>0 && i < length(T) %<= size(T,1) RBZ 271119 Changed to not route water past the last cell
    max_water=K_delta(i).*fieldCapacity(i);  %maximum amount of water (in m) that a grid cell can hold
    min_water = 0;  %minimum amount of water which stays in a cell (independent of soil type, but should be if "freezing = drying")
    
    actual_water= max( min_water,  wc(i).*K_delta(i)+dwc(i) );
    

    dwc(i+1)=dwc(i+1) + max(0, actual_water-max_water);  %when excess water, move it to next grid cell
    wc(i)=min(max_water, actual_water)./K_delta(i);  % in fraction
    i=i+1;
end

excess_water=dwc(i); %in m

i=i-1;

while i>=1 && excess_water>0
    max_water=K_delta(i).*porosity(i);
    actual_water=wc(i).*K_delta(i)+excess_water;
    wc(i)=min(actual_water, max_water)./K_delta(i);
    excess_water=max(0, actual_water-wc(i).*K_delta(i));
    i=i-1;
end


ground.STATVAR.water = wc .* ground.STATVAR.layerThick;
ground.STATVAR.waterIce = double(ground.STATVAR.T>0) .* ground.STATVAR.water + double(ground.STATVAR.T<=0) .* ground.STATVAR.waterIce;
ground.STATVAR.surface_runoff = double(excess_water>0).*excess_water; % surface runoff  if excess_water>0, could go to external reservoir