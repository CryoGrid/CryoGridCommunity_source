% inhertits from GROUND_freeW_bucketW and adds SEB as upper boundary
% no interaction with snow is possible here

classdef GROUND_freeW_bucketW_seb_snow < GROUND_freeW_bucketW_seb
    properties
        IA_CHILD
    end
    methods
        
        %mandatory functions for each class
        
        function xls_out = write_excel(ground)
            xls_out = {'CLASS','index',NaN,NaN,NaN;'GROUND_freeW_bucketW_seb_snow',1,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN;NaN,'value','default','unit',NaN;'albedo',0.200000000000000,0.150000000000000,'[-]','surface albedo';'epsilon',0.990000000000000,0.990000000000000,'[-]','surface emissivity';'z0',0.00100000000000000,0.00100000000000000,'[m]','roughness length';'rootDepth',0.200000000000000,NaN,'[m]',NaN;'evaporationDepth',0.0500000000000000,NaN,'[m]',NaN;'ratioET',0.500000000000000,NaN,'[-]',NaN;'hydraulicConductivity',1.00000000000000e-05,NaN,'[m/sec]',NaN;'    ','    ','    ',NaN,'    ';'dt_max',3600,3600,'[sec]','longest possible timestep';'dE_max',50000,50000,'[J/m3]','maximum change of energy per timestep';'CLASS_END',NaN,NaN,NaN,NaN};
        end

        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            ground = provide_variables@GROUND_freeW_bucketW_seb(ground); %call function of the base class
        end
        
        function ground = assign_global_variables(ground, forcing)
            ground = assign_global_variables@GROUND_freeW_bucketW_seb(ground, forcing); %call function of the base class
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            ground = initialize_STATVAR_from_file@GROUND_freeW_bucketW_seb(ground, grid, forcing, depths);
        end
        
        
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            ground = get_boundary_condition_u@GROUND_freeW_bucketW_seb(ground, forcing);       
            forcing.TEMP.rainfall = 0;
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = get_boundary_condition_u(ground.IA_CHILD, forcing); %call boundary condition for child
            end
                  
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
            ground = get_boundary_condition_l@GROUND_freeW_bucketW(ground, forcing);
        end
        
        
        function ground = get_derivatives_prognostic(ground)
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = get_derivative_energy(ground.IA_CHILD); % boundary condition for child, omacts ground.TEMP.F_ub
                
            else
                ground = get_derivatives_prognostic@GROUND_freeW_bucketW(ground);
            end
            
            ground = get_derivative_water(ground); 
            %add the water flux bc
           
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
           timestep = get_timestep@GROUND_freeW_bucketW(ground);
        end
        
        function ground = advance_prognostic(ground, timestep) 
            ground = advance_prognostic@GROUND_freeW_bucketW_seb(ground, timestep); %advance energy and route down water
            
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = advance_prognostic(ground.IA_CHILD, timestep); %call function for child CHECK!!
            end
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            ground = compute_diagnostic_first_cell@GROUND_freeW_bucketW_seb(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, forcing)
            ground = get_T_water(ground);
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = compute_diagnostic(ground.IA_CHILD); %call function for child
                ground.IA_CHILD = check_trigger(ground.IA_CHILD);
            end
            ground = conductivity(ground);
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = mix_conductivity(ground.IA_CHILD); %call function for child
            end
            
            %ground = compute_diagnostic@GROUND_freeW_bucketW(ground, forcing);
        end
        
        function ground = troubleshoot(ground)
            ground = checkNaN(ground);
        end
        
        %non-mandatory functions
        function ground = get_T_water(ground)
            ground = get_T_water@GROUND_base_class(ground);
        end
        
        function ground = conductivity(ground)
            ground = conductivity@GROUND_base_class(ground);
        end

        
        function ground = get_derivative_water(ground)
            ground = get_derivative_water@GROUND_freeW_bucketW(ground);
        end
    
    end
    
end
