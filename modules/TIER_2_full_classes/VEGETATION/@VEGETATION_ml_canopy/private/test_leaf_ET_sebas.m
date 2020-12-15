leaf_ET=[];
transp_flux_per_layer = [];

for i=1:size(out.FORCING,2)
    frac_sun_shade = cat(3, out.VEGETATION{1,i}.STATVAR.vegetation.flux.fracsun, out.VEGETATION{1,i}.STATVAR.vegetation.flux.fracsha);
    let=out.VEGETATION{1,i}.STATVAR.vegetation.mlcanopyinst.trleaf .* frac_sun_shade;
    let=sum(let,3); %per m2 leaf per m2 ground
    let= let .* out.VEGETATION{1,i}.STATVAR.vegetation.canopy.dlai; %per m2 ground
    let_tot = sum(0.018e-3 .* let); %in m3/me water/sec
    tfpl = let_tot .* out.VEGETATION{1,i}.STATVAR.vegetation.mlcanopyinst.soil_et_loss;
    transp_flux_per_layer = [transp_flux_per_layer; tfpl];
    leaf_ET=[leaf_ET; let];
end