% classdef CLM_Vegetation
%     properties
%         mlcanopyinst %constants 
%     end
%     
%     methods

        %mandatory functions for each module
%         vegetation = initialize_mlcanopyinst(FORCING);

        
%         %upper boundry condition
%         function forest = get_boundary_condition_u(forest, forcing) 
%             forest = surface_energy_balance(forest, forcing);
%             forest.TEMP.rainfall = forcing.TEMP.rainfall ./1000 ./(24.*3600);
%         end
%         
%         %lower boundry condition
%         function forest = get_boundary_condition_l(forest) 
%             forest.TEMP.F_lb = get_F_lb(forest);
%         end
%         
%         
%         function forest = get_derivatives_prognostic(forest)
%             forest.TEMP.d_energy = get_derivative_energy(forest);
%         end
%         
%         function timestep = get_timestep(forest)  %could involve check for several state variables
%              timestep = forest.PARA.d_e_max ./ (max(abs(forest.TEMP.d_energy) ./ forest.STATVAR.D));  
%              
%         end
%         
%         function forest = advance_prognostic(forest, timestep) %real timestep derived as minimum of several classes in [sec] here!
%             forest.STATVAR.energy = forest.STATVAR.energy + timestep .* forest.TEMP.d_energy;
%         end
%         
%         function forest = compute_diagnostic_first_cell(forest, forcing)
%             forest = L_star(forest, forcing);
%         end
%         
%         function forest = compute_diagnostic(forest)
%             forest = get_T_water(forest); 
%             forest.STATVAR.cond = conductivity_veg(forest);
% %             forest = L_star(forest, forcing);
%         end
% 
% %         function [trigger, forest] = check_triggers (forest, trigger)
%             if forest.STATVAR.swe > forest.PARA.swe_per_cell./2
%                 snow_energy = -forest.STATVAR.swe .* forest.CONST.Lf + forest.STATVAR.T(1) .* forest.STATVAR.swe .* forest.CONST.c_i;
%                 trigger=[trigger; {CREATE_SNOW(forest.STATVAR.swe, snow_energy, forest.STATVAR.T(1), forest.STATVAR.water_reservoir, forest.STATVAR.Lstar)}];
%                 forest.STATVAR.energy(1) = forest.STATVAR.energy(1) - snow_energy;
%                 forest.STATVAR.swe = 0;
%                 forest.STATVAR.water_reservoir = 0; 
%                 forest.TEMP.snowfall=0;
%                 forest.TEMP.rainfall = 0;
%                 forest.TEMP.snow_energy = 0;
%             end
  