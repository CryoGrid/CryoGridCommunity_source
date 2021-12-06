function TC = ThermalConductivity(WaterIce, Water, Mineral, Organic)

ka = 0.025;      %air [Hillel(1982)]
kw = 0.57;        %water [Hillel(1982)]
ko = 0.25;        %organic [Hillel(1982)]
km = 3.8;         %mineral [Hillel(1982)]
%km = 2.8;         %average Granite [Bejan and Kraus (2003) in Dong et al. (2015)]
ki = 2.2;         %ice [Hillel(1982)]

n = length(WaterIce);
TC = zeros(n,1);
for i=1:n
    ice = WaterIce(i) - Water(i);
    air = 1.0 - WaterIce(i) - Mineral(i) - Organic(i);
    TC(i,1) = (Water(i).* kw.^0.5 + ice.* ki.^0.5 + Mineral(i).* km.^0.5 + Organic(i).* ko.^0.5 + air.* ka.^0.5).^2.0;
end

end
