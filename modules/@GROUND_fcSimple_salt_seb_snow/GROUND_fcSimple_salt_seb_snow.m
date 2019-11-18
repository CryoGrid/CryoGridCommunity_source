% inherits from GROUND_fcSimple_salt_seb and adds interaction with a SNOW
% class through the IA_SNOW_GROUND_fcSimple class

classdef GROUND_fcSimple_salt_seb_snow < GROUND_fcSimple_salt_seb
    properties
        IA_CHILD
    end
    
    
    methods

        %mandatory functions for each class
        
        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            ground = provide_variables@GROUND_fcSimple_salt_seb(ground);
        end
        
        function variable = initialize_from_file(ground, variable, section)
            variable = initialize_from_file@GROUND_fcSimple_salt_seb(ground, variable, section);
        end
        
        function ground = assign_global_variables(ground, forcing)
            ground = assign_global_variables@GROUND_fcSimple_salt_seb(ground, forcing);
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            ground = initialize_STATVAR_from_file@GROUND_fcSimple_salt_seb(ground, grid, forcing, depths);
        end
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            ground = get_boundary_condition_u@GROUND_fcSimple_salt_seb(ground, forcing);
            
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = get_boundary_condition_u(ground.IA_CHILD, forcing); %call boundary condition for child
            end
        end
        
        function ground = get_boundary_condition_l(ground) 
            ground = get_boundary_condition_l@GROUND_fcSimple_salt_seb(ground);
        end
        
        
        function ground = get_derivatives_prognostic(ground)
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = get_derivative_energy(ground.IA_CHILD); % boundary condition for child, omacts ground.TEMP.F_ub
            else
                ground = get_derivative_energy(ground);
            end
            ground = get_derivative_salt(ground); 
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables, in this case only check energy derivative
               timestep = get_timestep@GROUND_fcSimple_salt_seb(ground);
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            ground = advance_prognostic@GROUND_fcSimple_salt_seb(ground, timestep);
            
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = advance_prognostic(ground.IA_CHILD, timestep); %call function for child
            end
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing);
            ground = compute_diagnostic_first_cell@GROUND_fcSimple_salt_seb(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, forcing)

            ground = get_T_water_salt_FreezeDepress_Xice(ground);
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = compute_diagnostic(ground.IA_CHILD); %call function for child   -> modify this function so that SWE is give out
                ground.IA_CHILD = check_trigger(ground.IA_CHILD);
            end
            ground = conductivity(ground);
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = mix_conductivity(ground.IA_CHILD); %call function for child
            end
            
            ground = diffusivity_salt(ground); 
        end
        
        
        
        %non-mandatory functions -> required here so that they are usable
        %in subclasses

        function ground = get_T_water_salt_FreezeDepress_Xice(ground)
            ground = get_T_water_salt_FreezeDepress_Xice@GROUND_fcSimple_salt(ground);
        end
        
        function ground = get_derivative_salt(ground)
            ground = get_derivative_salt@GROUND_fcSimple_salt(ground);
        end

        function ground = conductivity(ground)
            ground = conductivity@GROUND_base_class(ground);
        end
        
        function  ground = get_derivative_energy(ground)
            ground = get_derivative_energy@GROUND_base_class(ground);
        end
        
        
    end
end