%========================================================================
% CryoGrid TIER1 library class for functions related heat conduction
% using only an upper boundary temperature forcing.
%
% T. Ingeman-Nielsen, S. Westermann, December 2021
%========================================================================

classdef LB_TEMPERATURE_FORCING < BASE
    
    methods
        %------boundary conditions-----
        function ground = get_boundary_condition_l(ground, forcing)
            % Calculate upper boundary flux and energy change due to
            % temperature forcing.
            ground.TEMP.F_lb = -(ground.STATVAR.T(end,end)-forcing.TEMP.T_lb).*ground.STATVAR.thermCond(end,end)./(ground.STATVAR.layerThick(end,end)/2); 
            ground.TEMP.d_energy(end,end) = ground.TEMP.d_energy(end,end) + ground.TEMP.F_lb.*ground.STATVAR.area(end,end);
        end        
        
    end
end

