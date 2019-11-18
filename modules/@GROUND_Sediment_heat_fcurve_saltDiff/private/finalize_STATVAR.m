function ground = finalize_STATVAR(ground)
    %porosity, Tmelt, initial conditions are already claculated in the base
    %class
    
    %conductivity, heat capacity and liquid water content
    ground = getThermalProps_wSalt(ground);
end
