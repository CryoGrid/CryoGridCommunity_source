% inhertits from GROUND_base_class and adds a bucket scheme for water

classdef GROUND_freeW_bucketW < GROUND_base_class

    
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
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            ground = initialize_STATVAR_from_file@GROUND_base_class(ground, grid, forcing, depths);
            ground = finalize_STATVAR(ground); %assign all variables, that must be calculated or assigned otherwise
        end
        

        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            ground.TEMP.F_ub_water = forcing.TEMP.rainfall ./ 1000 ./ 24 ./3600;  %possibly add water from external source here 
            ground.TEMP.T_rainWater =  max(0,forcing.TEMP.Tair);
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
            ground = get_boundary_condition_l@GROUND_base_class(ground, forcing);
            ground.TEMP.F_lb_water = 0; % zero flux if used as bottom class
        end
        
        
        function ground = get_derivatives_prognostic(ground)
            ground = get_derivatives_prognostic@GROUND_base_class(ground);
            ground = get_derivative_water(ground); %NEW
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
           timestep = get_timestep@GROUND_base_class(ground);
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            ground = advance_prognostic@GROUND_base_class(ground, timestep);
            ground = advance_prognostic_water(ground, timestep); %NEW
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            %empty here 
        end
        
        function ground = compute_diagnostic(ground, forcing)
            ground = compute_diagnostic@GROUND_base_class(ground, forcing);
        end
        
        
        %non-mandatory functions
        function ground = get_derivative_water(ground)
          
            saturation = ground.STATVAR.water ./ (ground.STATVAR.layerThick - ground.STATVAR.mineral - ground.STATVAR.organic);
            waterMobile = double(saturation > ground.STATVAR.field_capacity);
            ground.TEMP.d_water_out = waterMobile .* ground.PARA.hydraulicConductivity .* ground.STATVAR.water ./ ground.STATVAR.layerThick;
            ground.TEMP.d_water_in = ground.TEMP.d_water_out .*0;
            ground.TEMP.d_water_in(2:end,1) = ground.TEMP.d_water_out(1:end-1,1);
            ground.TEMP.d_water_in(1,1) = ground.TEMP.F_ub_water;
            if ~isempty(ground.IA_NEXT)
                get_boundary_condition_water_m(ground.IA_NEXT);
            else
                ground.TEMP.d_water_out(end,1) = ground.TEMP.F_lb_water;
            end
        end
        
        function ground = advance_prognostic_water(ground, timestep)
            %ground.STATVAR.water  = ground.STAVAR.water + ground.TEMP.d_water_ET .* timestep; %add this function when SEB is present

            
            ground.TEMP.d_water_in = ground.TEMP.d_water_in .* timestep;
            ground.TEMP.d_water_out = ground.TEMP.d_water_out .* timestep;
            %limit outflow to field capacity
            ground.TEMP.d_water_out  = min(ground.TEMP.d_water_out, max(0, ground.STATVAR.water - ground.STATVAR.field_capacity .* (ground.STATVAR.layerThick - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.ice)));
            ground.TEMP.d_water_in(2:end,1) = ground.TEMP.d_water_out(1:end-1,1);
            %limit inflow so that unity is not exceeded
            ground.TEMP.d_water_in = min(ground.TEMP.d_water_in, ground.STATVAR.layerThick - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce);
            ground.TEMP.d_water_out(1:end-1,1) = ground.TEMP.d_water_in(2:end,1);
            
            finalize_boundary_condition_water_m(ground.IA_NEXT, timestep);
            
            energy_out = ground.CONST.c_w .* ground.STATVAR.T .* ground.TEMP.d_water_out;
            energy_in = energy_out.*0;
            energy_in (2:end,1) = energy_out(1:end-1,1);
            energy_in(1,1) = ground.CONST.c_w .* ground.TEMP.T_rainWater .* ground.TEMP.d_water_in(1,1);
            
            ground.STATVAR.waterIce = ground.STATVAR.waterIce - ground.TEMP.d_water_out + ground.TEMP.d_water_in;
            ground.STATVAR.energy = ground.STATVAR.energy -energy_out + energy_in;

        end
        
        
    end
    
end
