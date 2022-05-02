%========================================================================
% CryoGrid TIER1 library class for functions related heat conduction
% using only an upper boundary temperature forcing.
%
% T. Ingeman-Nielsen, S. Westermann, December 2021
%========================================================================

classdef UB_TEMPERATURE_FORCING < BASE
    
    methods
        %------boundary conditions-----
        function ground = get_boundary_condition_u(ground, forcing)
            % Calculate upper boundary flux and energy change due to
            % temperature forcing.
            ground.TEMP.F_ub = -(forcing.TEMP.Tair - ground.STATVAR.T(1,1)).*ground.STATVAR.thermCond(1,1)./(ground.STATVAR.layerThick(1,1)/2); 
            ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) - ground.TEMP.F_ub.*ground.STATVAR.area(1,1);
        end        
        
        function ground = get_ub_temperature(ground, forcing)
            % Calculate upper boundary flux and energy change due to
            % temperature forcing.
            ground.TEMP.F_ub = -(forcing.TEMP.Tair - ground.STATVAR.T(1,1)).*ground.STATVAR.thermCond(1,1)./(ground.STATVAR.layerThick(1,1)/2); 
            ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) - ground.TEMP.F_ub.*ground.STATVAR.area(1,1);
        end  
        
    end
end

