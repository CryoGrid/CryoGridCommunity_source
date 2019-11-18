% inhertits from GROUND_base_class and adds SEB as upper boundary
% no interaction with snow is possible here

classdef GROUND_freeW_seb < GROUND_base_class
    
    methods
        
        %mandatory functions for each class

        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            
            ground = provide_variables@GROUND_base_class(ground); %call function of the base class
            ground = provide_PARA(ground); %add additional variables
            ground = provide_CONST(ground);
            ground = provide_STATVAR(ground);
        end
        
        function ground = assign_global_variables(ground, forcing)
            ground = assign_global_variables@GROUND_base_class(ground, forcing); %call function of the base class
            ground.PARA.airT_height = forcing.PARA.airT_height;
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            ground = initialize_STATVAR_from_file@GROUND_base_class(ground, grid, forcing, depths);

            ground = finalize_STATVAR(ground); %assign all variables, that must be calculated or assigned otherwise
        end
        
               
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            ground = surface_energy_balance(ground, forcing);      
            
%             figure(20)
%             plot(ground.STATVAR.T)
%             title('Ground temp')
%             xlabel('Ground layer')
%             ylabel('Temperature')
%             hold on
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
            ground = get_boundary_condition_l@GROUND_base_class(ground, forcing);
            
            figure(20)
            plot(ground.STATVAR.T)
            title('Ground temp')
            xlabel('Ground layer')
            ylabel('Temperature')
            hold on
        end
        
        
        function ground = get_derivatives_prognostic(ground)
            ground = get_derivatives_prognostic@GROUND_base_class(ground);
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
           timestep = get_timestep@GROUND_base_class(ground);
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            ground = advance_prognostic@GROUND_base_class(ground, timestep);
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, forcing)
            ground = compute_diagnostic@GROUND_base_class(ground, forcing);
        end

%         non-mandatory fucntions
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe = Q_eq(ground, forcing);
            
            ground.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;

        end

    end
    
end
