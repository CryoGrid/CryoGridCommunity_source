% inhertits from GROUND_fcSimple_salt and adds SEB as upper boundary
% no interaction with snow is possible here

classdef GROUND_fcSimple_salt_seb < GROUND_fcSimple_salt
   
    
    methods

        %mandatory functions for each class
        
        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            ground = provide_variables@GROUND_fcSimple_salt(ground);
            ground = provide_PARA(ground);
            ground = provide_CONST(ground);
            ground = provide_STATVAR(ground);
        end
        
        function variable = initialize_from_file(ground, variable, section)
            variable = initialize_from_file@GROUND_fcSimple_salt(ground, variable, section);
        end
        
        function ground = assign_global_variables(ground, forcing)
            ground = assign_global_variables@GROUND_fcSimple_salt(ground, forcing);
            ground.PARA.airT_height = forcing.PARA.airT_height;
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            ground = initialize_STATVAR_from_file@GROUND_fcSimple_salt(ground, grid, forcing, depths);
            ground = finalize_STATVAR(ground); %assign all variables, that must be calculated or assigned otherwise
        end
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            ground = surface_energy_balance(ground, forcing);
            ground.TEMP.saltConcFlux_ub=0; %zero flux bc, could be changed to Dirichlet BC, but this makes only sense when temperature bc is also Dirichlet             
        end
        
        function ground = get_boundary_condition_l(ground) 
            ground = get_boundary_condition_l@GROUND_fcSimple_salt(ground);
            ground.TEMP.saltConcFlux_lb=0; %zero flux bc (default, but could be changed)
        end
        
        function ground = get_derivatives_prognostic(ground)
            ground = get_derivatives_prognostic@GROUND_fcSimple_salt(ground);
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables, in this case only check energy derivative
               timestep = get_timestep@GROUND_fcSimple_salt(ground);
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            ground = advance_prognostic@GROUND_fcSimple_salt(ground, timestep);
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing);
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, forcing)
            ground = compute_diagnostic@GROUND_fcSimple_salt(ground, forcing);
        end
        

        %non-mandatory functions -> required here so that they are usable
        %in subclasses
                
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe = Q_eq(ground, forcing);
            ground.TEMP.heatFlux_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
        end

        
    end
end