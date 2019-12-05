    function conductivity2 = conductivity2(water, ice, mineral, organic, ground)
        
        ka = ground.CONST.k_a; %0.025;       %air [Hillel(1982)]
        kw = ground.CONST.k_w; %0.57;        %water [Hillel(1982)]
        ko = ground.CONST.k_o; %0.25;        %organic [Hillel(1982)]
        km = ground.CONST.k_m; %soil.kh_bedrock;     %mineral
        ki = ground.CONST.k_i; %2.2;         %ice [Hillel(1982)]
        
        air=1-water-ice-mineral-organic;
        
        conductivity2= (water.* kw.^0.5 + ice.* ki.^0.5 + mineral.* km.^0.5 + organic.* ko.^0.5 + air.* ka.^0.5).^2;
    end