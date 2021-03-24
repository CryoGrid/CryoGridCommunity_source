function [phi_m] = phi_m_monin_obukhov (zeta)

% --- Evaluate the Monin-Obukhov phi function for momentum at x

if (zeta < 0)
   phi_m = (1 - 16 * zeta)^(-0.25);
else
   phi_m = 1 + 5 * zeta;
end
