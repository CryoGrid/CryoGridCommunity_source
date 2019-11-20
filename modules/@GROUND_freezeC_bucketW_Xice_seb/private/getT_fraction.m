function fraction = getT_fraction(ground)

T = ground.STATVAR.T;
wc = ground.STATVAR.water ./ ground.STATVAR.layerThick;
start_reduction = ground.STATVAR.field_capacity;

fraction=double(T>0).*(double(wc>=start_reduction) + double(wc<start_reduction).*0.25.*(1-cos(pi().*wc./start_reduction)).^2);
