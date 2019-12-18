function [ lateral, waterflux ] = lateral_water_oscillations( lateral, waterflux)
% Check if the lateral fluxes generate oscillations in the water content.
% So identify alternating patterns and brakes it.

pattern=lateral.TEMP.lastWaterChange;
pattern=circshift(pattern,[0,-1]);
pattern(end)=waterflux;

if sum(isnan(pattern))==0
    pos=sum(find(pattern>0));
    if mod(pos,2)==0 && pos>0 && pos<6
        waterflux=waterflux/3; % 3 seems to perform well, below, some oscillations are not killed, above, oscillations like --+--+--+--+--+ start appearing
        pattern(end)=waterflux;
        disp('lateral_water_oscillations.m : dividing water flux')
    end
end

lateral.TEMP.lastWaterChange=pattern;

end