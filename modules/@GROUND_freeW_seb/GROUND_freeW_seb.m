% inhertits from GROUND_base_class and adds SEB as upper boundary
% no interaction with snow is possible here

classdef GROUND_freeW_seb < GROUND_base_class

    
    
    methods
        
        %mandatory functions for each class
        
        function xls_out = write_excel(ground) %xls_out is a cell array corresponding to the class-specific content of the parameter excel file
            xls_out = {'CLASS','index',NaN,NaN,NaN;'GROUND_freeW_seb',1,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN;NaN,'value','default','unit',NaN;'albedo',0.150000000000000,0.150000000000000,'[-]','surface albedo';'epsilon',0.990000000000000,0.990000000000000,'[-]','surface emissivity';'rs',100,100,'[-]','surface resistance against evapotransipration';'z0',0.00100000000000000,0.00100000000000000,'[m]','roughness length';'    ','    ','    ',NaN,'    ';'dt_max',3600,3600,'[sec]','longest possible timestep';'dE_max',50000,50000,'[J/m3]','maximum change of energy per timestep';'CLASS_END',NaN,NaN,NaN,NaN};
        end

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
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
            ground = get_boundary_condition_l@GROUND_base_class(ground, forcing);
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
        
        
        %non-mandatory fucntions
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe = Q_eq(ground, forcing);
            
            ground.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
        end
    
    end
    
end
