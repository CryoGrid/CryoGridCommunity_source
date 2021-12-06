function k_eff_ground = k_eff_compiled(k, K_delta)

k_eff_ground = k(4:end-1,:).*k(5:end,:).*(K_delta(5:end-1,:)./2 + K_delta(6:end,:)./2) ./ (k(4:end-1,:).*K_delta(6:end,:)./2 + k(5:end,:).*K_delta(5:end-1,:)./2 ); %size N

