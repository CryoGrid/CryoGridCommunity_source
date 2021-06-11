% function res   = dE_dt_compiled(k_eff, T, cT_delta)
% 
%          res = (k_eff(2:end,:).*(T(3:end,:)-T(2:end-1,:))./cT_delta(2:end,:) -...
%                 k_eff(1:end-1,:).*(T(2:end-1,:)-T(1:end-2,:))./cT_delta(1:end-1,:));  %size N-1
           
            
function d_energy   = dE_dt_compiled(d_energy, thermCond_eff, T, layerDistance)

%downwards flux
d_energy = d_energy - thermCond_eff.*(T(2:end,:) - T(1:end-1,:)) ./ layerDistance;
%upwards flux, lower boundary already added
d_energy(1:end-1,:) = d_energy(1:end-1,:) + thermCond_eff(2:end,:).*(T(3:end,:) - T(2:end-1,:)) ./ layerDistance(2:end,:);

