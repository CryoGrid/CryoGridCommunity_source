function  thermCond_eff = thermCond_eff_compiled(thermCond, layerThick)

thermCond_eff = thermCond(5:end-1,:).*thermCond(6:end,:) .* (layerThick(5:end-1,:)./2 + layerThick(6:end,:)./2) ./ (thermCond(5:end-1,:).*layerThick(6:end,:)./2 + thermCond(6:end,:).*layerThick(5:end-1,:)./2 ); %size N
