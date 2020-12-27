function thermCond = thermalConductivity_CLM4_5(mineral, organic, waterIce, water, ice, T)

k_m = 3; %can be expressed in terms of sand 8.8 and clay contents 2.9W/mK
k_i = 2.2;
k_w = 0.57;
k_o = 0.25;
k_dry_organic = 0.05; %slightly nonsense...





porosity = 1 - mineral - organic;
organic_fraction = organic ./ (mineral + organic);

k_solids = organic_fraction .* k_o + (1- organic_fraction) .* k_m;
k_sat = k_solids.^(1-porosity) .* k_w .^(water./waterIce.* porosity) .* k_i .^(ice./waterIce.* porosity); 

bulk_density = 2700 .* (1-porosity);
k_dry_mineral = (0.135 .* bulk_density + 64.7) ./ (2700 - 0.947 .* bulk_density);
k_dry = organic_fraction .* k_dry_organic + (1- organic_fraction) .* k_dry_mineral;

saturation = waterIce./porosity;
Kersten_number = double(T>=0) .* max(0, log(saturation) ./ log(10) + 1) +  double(T<0) .* saturation;

thermCond = Kersten_number .* k_sat + (1- Kersten_number) .* k_dry;



