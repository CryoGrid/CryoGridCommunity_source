function  thermCond = thermCond_compiled(T, T_end_freezing, k_frozen, k_freezing, k_thawed)

thermCond = double(T(5:end,:) < T_end_freezing) .* k_frozen + double(T(5:end,:) > 0) .* k_thawed + ...
    double(T(5:end,:) >= T_end_freezing & T(5:end,:) <= 0) .* k_freezing;
                