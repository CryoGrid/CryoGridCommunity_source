%========================================================================
% CryoGrid TIER1 library class for functions related to air convection
% CAUTION: this is highly experimental and should not be used!!
% The theory will likely not withstand a peer-review process!
% S. Westermann, October 2020
%========================================================================


classdef AIR_CONVECTION < BASE
    
    methods
        
        function ground = get_boundary_condition_u_convection(ground, forcing)
            
            rho_air = ground.PARA.pressure ./ (ground.CONST.R_spec.*([forcing.TEMP.Tair; ground.STATVAR.T(1,1)] + ground.CONST.Tmfw));
            
            delta_p = double(rho_air(1,1) > rho_air(2,1)) .* (rho_air(1,1) - rho_air(2,1)) ./ 4 .* ground.CONST.g ;

            %rough pipe turbulent regime, constant Darcy friction factor
            velocity_turbulent = sqrt(2 .* ground.STATVAR.diamater_pipe(1,1) ./ mean(rho_air,1) ./ground.CONST.Darcy_friction_factor .* delta_p ./ (ground.STATVAR.layerThick(1,1)./2)); %Darcy Weisbach   

            %laminar flow
            velocity_laminar = ground.STATVAR.diamater_pipe(1,1).^2 ./ (32.*ground.CONST.viscosity_air) .* delta_p ./ (ground.STATVAR.layerThick(1,1)./2);
           
            air_flux = pi() ./ 4 .* ground.STATVAR.diamater_pipe(1,1).^2 .* min(velocity_turbulent, velocity_laminar)  .* ground.STATVAR.number_of_pipes(1,1);        
     
            energy_convection_down =  air_flux.* ground.CONST.cp .*(rho_air(1,1) + rho_air(2,1))./2 .* forcing.TEMP.Tair .* ground.STATVAR.area(1,1);
            energy_convection_up = air_flux.* ground.CONST.cp .* (rho_air(1,1) + rho_air(2,1))./2 .* ground.STATVAR.T(1,1) .* ground.STATVAR.area(1,1);
            
            ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) - energy_convection_up + energy_convection_down;
        end
        
        %-----derivatives----------
        
        function ground = get_derivative_air_convection_Darcy(ground) %dont use, this assumes laminar flow and overestimates flow at high block sizes

            permeability_air = ground.STATVAR.permeability_air(1:end-1,1) .* ground.STATVAR.permeability_air(2:end,1) ./ ...
                (ground.STATVAR.permeability_air(1:end-1,1) .* ground.STATVAR.layerThick(2:end,1) ./2 + ground.STATVAR.permeability_air(2:end,1) .* ground.STATVAR.layerThick(1:end-1,1) ./ 2);

            permeability_air(isnan(permeability_air)) = 0;
            
            rho_air = ground.PARA.pressure ./ (ground.CONST.R_spec.*(ground.STATVAR.T + ground.CONST.Tmfw));
            
            delta_p = double(rho_air(1:end-1,1) > rho_air(2:end,1)) .* (rho_air(1:end-1,1) - (rho_air(1:end-1,1)+rho_air(2:end,1))./2) .* ground.STATVAR.layerThick(2:end) ./ 2 .* ground.CONST.g ;

            air_flux =  permeability_air ./ ground.CONST.viscosity_air .* delta_p; %m3/sec m2
                        
            energy_convection_down =  air_flux.* ground.CONST.cp .*(rho_air(1:end-1,1)+rho_air(2:end,1))./2 .* ground.STATVAR.T(1:end-1,1) .* ground.STATVAR.area(1:end-1,1);
            energy_convection_up = air_flux.* ground.CONST.cp .* (rho_air(1:end-1,1)+rho_air(2:end,1))./2 .* ground.STATVAR.T(2:end,1) .* ground.STATVAR.area(2:end,1);
            
            ground.TEMP.d_energy(1:end-1,1) = ground.TEMP.d_energy(1:end-1,1) + energy_convection_up - energy_convection_down;
            ground.TEMP.d_energy(2:end,1) = ground.TEMP.d_energy(2:end,1) - energy_convection_up + energy_convection_down;
        end
        
                
        function ground = get_derivative_air_convection_Darcy_Weisbach(ground)
            
            rho_air = ground.PARA.pressure ./ (ground.CONST.R_spec.*(ground.STATVAR.T + ground.CONST.Tmfw));
            
            delta_p = double(rho_air(1:end-1,1) > rho_air(2:end,1)) .* (rho_air(1:end-1,1) - rho_air(2:end,1)) .* ground.STATVAR.layerThick(1:end-1,1).* ground.STATVAR.layerThick(2:end,1) ./ (ground.STATVAR.layerThick(1:end-1,1) + ground.STATVAR.layerThick(2:end,1)).^2 .* ground.CONST.g; 

            %rough pipe turbulent regime, constant Darcy friction factor
            delta_p_upperCell = ground.STATVAR.diamater_pipe(2:end,1).^5 .* ground.STATVAR.number_of_pipes(2:end,1).^2 .* delta_p ./ ground.STATVAR.layerThick(2:end,1) .* ...
                (ground.STATVAR.diamater_pipe(2:end,1).^5 .* ground.STATVAR.number_of_pipes(2:end,1).^2 ./ ground.STATVAR.layerThick(2:end,1) +  ground.STATVAR.diamater_pipe(1:end-1,1).^5 .* ground.STATVAR.number_of_pipes(1:end-1,1).^2 ./ ground.STATVAR.layerThick(1:end-1,1)).^-1;
            delta_p_upperCell(isnan(delta_p_upperCell)) = 0;
            
            velocity_upper_cell_turbulent = sqrt(2 .* ground.STATVAR.diamater_pipe(1:end-1,1) ./ mean(rho_air,1) ./ground.CONST.Darcy_friction_factor .* delta_p_upperCell ./ (ground.STATVAR.layerThick(1:end-1,1)./2));%Darcy Weisbach   

            %laminar flow
            delta_p_upperCell = ground.STATVAR.diamater_pipe(2:end,1).^4 .* ground.STATVAR.number_of_pipes(2:end,1) .* delta_p ./ ground.STATVAR.layerThick(2:end,1) .* ...
                (ground.STATVAR.diamater_pipe(2:end,1).^4 .* ground.STATVAR.number_of_pipes(2:end,1) ./ ground.STATVAR.layerThick(2:end,1) +  ground.STATVAR.diamater_pipe(1:end-1,1).^4 .* ground.STATVAR.number_of_pipes(1:end-1,1) ./ ground.STATVAR.layerThick(1:end-1,1)).^-1;
            delta_p_upperCell(isnan(delta_p_upperCell)) = 0;
            
            velocity_upper_cell_laminar = ground.STATVAR.diamater_pipe(1:end-1,1).^2 ./ (32.*ground.CONST.viscosity_air) .* delta_p_upperCell ./ (ground.STATVAR.layerThick(1:end-1,1)./2);
            
            
            ground.TEMP.air_flux = pi() ./ 4 .* ground.STATVAR.diamater_pipe(1:end-1).^2 .* min(velocity_upper_cell_turbulent, velocity_upper_cell_laminar)  .* ground.STATVAR.number_of_pipes(1:end-1);        
            
            energy_convection_down =  ground.TEMP.air_flux.* ground.CONST.cp .*(rho_air(1:end-1,1)+rho_air(2:end,1))./2 .* ground.STATVAR.T(1:end-1,1) .* ground.STATVAR.area(1:end-1,1);
            energy_convection_up = ground.TEMP.air_flux.* ground.CONST.cp .* (rho_air(1:end-1,1)+rho_air(2:end,1))./2 .* ground.STATVAR.T(2:end,1) .* ground.STATVAR.area(2:end,1);
            
            ground.TEMP.d_energy(1:end-1,1) = ground.TEMP.d_energy(1:end-1,1) + energy_convection_up - energy_convection_down;
            ground.TEMP.d_energy(2:end,1) = ground.TEMP.d_energy(2:end,1) - energy_convection_up + energy_convection_down;
        end
        
        
        function timestep = get_timestep_air_convection(ground) %leads to super-small timesteps
            air_flux_nonzero = ground.TEMP.air_flux > 0;
            timestep_up = min((ground.STATVAR.layerThick(air_flux_nonzero) ./ ground.STATVAR.area(air_flux_nonzero) - ground.STATVAR.waterIce(air_flux_nonzero) - ground.STATVAR.mineral(air_flux_nonzero) - ground.STATVAR.organic(air_flux_nonzero)) ./ 2 ./ ground.TEMP.air_flux(air_flux_nonzero)); %no more than half the void space up

            timestep_down  = (ground.STATVAR.layerThick(1:end-1,1) ./ ground.STATVAR.area(1:end-1,1) - ground.STATVAR.waterIce(1:end-1,1) - ground.STATVAR.mineral(1:end-1,1) - ground.STATVAR.organic(1:end-1,1)) ./ 2 ./ ground.TEMP.air_flux(2:end,1);
            
            timestep_down = nanmin(timestep_down);
            timestep = min(timestep_up, timestep_down);
            timestep = max(1, timestep);
        end
        
        
        %---permeability air--------------
        function ground = pipes_Darcy_Weisbach(ground)
            ground.STATVAR.porosity = 1 - (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick ./ ground.STATVAR.area; 
            
            ground.STATVAR.diamater_pipe = 2./3 .*ground.STATVAR.porosity ./(1-ground.STATVAR.porosity) .* ground.STATVAR.grain_size;
            ground.STATVAR.number_of_pipes = ground.STATVAR.porosity ./ (pi()./4 .* ground.STATVAR.diamater_pipe.^2 .* ground.CONST.tortuosity_air); 
            ground.STATVAR.number_of_pipes(isnan(ground.STATVAR.number_of_pipes)) = 0; %if diameter is zero
        end
        
        function ground = permeability_air_Carman_Kozeny(ground)
            porosity = 1 - (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ground.STATVAR.permeability_air =  ground.STATVAR.grain_size.^2 ./ 180 .* porosity .^3 ./ (1 - porosity).^2;
            ground.STATVAR.porosity = porosity;
            
        end
        
        function ground = permeability_air_Rumpf_Gupte(ground)
            porosity = 1 - (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ground.STATVAR.permeability_air = grain_size.^2 ./ 5.6 .* porosity .^5.5 ./ 1.05; %Rumpf-Gupte model [m2]
        end
        
    end
end

