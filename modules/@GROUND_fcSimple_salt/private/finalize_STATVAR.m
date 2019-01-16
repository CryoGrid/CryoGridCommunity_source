function ground = finalize_STATVAR(ground)


ground = get_E_water_salt_FreezeDepress_Xice(ground); %energy, water, ice, salt_c_brine

ground = conductivity(ground);

ground = diffusivity_salt(ground); % [m2/sec]

