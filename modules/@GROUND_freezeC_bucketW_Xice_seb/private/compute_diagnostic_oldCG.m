function ground = compute_diagnostic_oldCG(ground)

T = ground.STATVAR.T;

cT_frozen = ground.LOOKUP.cT_frozen;
cT_thawed = ground.LOOKUP.cT_thawed;
capacity = ground.LOOKUP.capacity;
K_frozen = ground.LOOKUP.cT_frozen; 
K_thawed = ground.LOOKUP.cT_thawed;
conductivity = ground.LOOKUP.conductivity; 
arraySizeT = ground.PARA.arraySizeT;
liquidWaterContent = ground.LOOKUP.liquidWaterContent;
                                               

a=(T-cT_frozen)./(cT_thawed-cT_frozen)*(arraySizeT-2)+1; %T and c information live on same grid
posT=round((a<=1).*(-a+1)+a+(a>arraySizeT-1).*(arraySizeT-a));
posT(posT==0)=1;
posT(isnan(posT))=arraySizeT;


r=[1:size(capacity,1)]';
c=posT;
indices=(size(capacity,1))*(c-1)+r;
c_temp=capacity(indices);
lwc_temp=liquidWaterContent(indices);

a=(T-K_frozen)./(K_thawed-K_frozen)*(arraySizeT-2)+1;
posT=round((a<=1).*(-a+1)+a+(a>arraySizeT-1).*(arraySizeT-a));
posT(posT==0)=1;
posT(isnan(posT))=arraySizeT;


r=[1:size(conductivity,1)]';
c=posT;
indices=(size(conductivity,1))*(c-1)+r;
k_eff=conductivity(indices);

ground.STATVAR.heatCapacity = c_temp;
ground.STATVAR.thermCond = k_eff;

ground.STATVAR.water = double(ground.STATVAR.T<=0) .* lwc_temp .* ground.STATVAR.layerThick + double(ground.STATVAR.T>0) .* ground.STATVAR.water;
ground.STATVAR.ice = ground.STATVAR.waterIce - ground.STATVAR.water;
