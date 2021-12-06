function HC = HeatCapacity(WaterIce, Water, Mineral, Organic)

cw = 4.2*10^6; %[J/m^3K] heat capacity water
co = 2.5*10^6; %[J/m^3K]  heat capacity organic
cm = 2.0*10^6; %[J/m^3K]  heat capacity mineral
ca = 0.00125*10^6;%[J/m^3K]  heat capacity pore space
ci = 1.9*10^6;%[J/m^3K]  heat capacity ice

n = length(WaterIce);
HC = zeros(n,1);
for i=1:n
    air = 1.0 - WaterIce(i) - Mineral(i) - Organic(i);
    ice = WaterIce(i) - Water(i);
    HC(i,1) = Water(i).*cw + ice.*ci + Mineral(i).*cm + Organic(i).*co + air.*ca;
end

end
