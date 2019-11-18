% inhertits from GROUND_freeW_bucketW and adds SEB as upper boundary
% no interaction with snow is possible here

classdef GROUND_freeW_bucketW_seb < GROUND_freeW_bucketW

    
    methods
        
        %mandatory functions for each class

        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            
            ground = provide_variables@GROUND_freeW_bucketW(ground); %call function of the base class
            ground = provide_PARA(ground); %add additional variables
            ground = provide_CONST(ground);
            ground = provide_STATVAR(ground);
        end
        
        function ground = assign_global_variables(ground, forcing)
            ground = assign_global_variables@GROUND_freeW_bucketW(ground, forcing); %call function of the base class
            ground.PARA.airT_height = forcing.PARA.airT_height;
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            ground = initialize_STATVAR_from_file@GROUND_freeW_bucketW(ground, grid, forcing, depths);
            ground = finalize_STATVAR(ground); %assign all variables, that must be calculated or assigned otherwise
        end
        
        
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            ground = get_boundary_condition_u@GROUND_freeW_bucketW(ground, forcing);
            ground = surface_energy_balance(ground, forcing);
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
            ground = get_boundary_condition_l@GROUND_freeW_bucketW(ground, forcing);
        end
        
        
        function ground = get_derivatives_prognostic(ground)
            ground = get_derivatives_prognostic@GROUND_freeW_bucketW(ground);
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
           timestep = get_timestep@GROUND_freeW_bucketW(ground);
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            ground.STATVAR.waterIce = ground.STATVAR.waterIce + ground.TEMP.d_water_ET .* timestep; %subtract water from ET
            ground.STATVAR.energy = ground.STATVAR.energy + ground.CONST.c_w .* ground.STATVAR.T .* ground.TEMP.d_water_ET .* timestep; %adjust energy
            ground = advance_prognostic@GROUND_freeW_bucketW(ground, timestep); %advance energy and route down water
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, forcing)
            ground = compute_diagnostic@GROUND_freeW_bucketW(ground, forcing);
        end
        
        
        %non-mandatory functions
        function ground = surface_energy_balance(ground, forcing)
            ground.STATVAR.Lout = (1-ground.PARA.epsilon) .* forcing.TEMP.Lin + ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+ 273.15).^4;
            ground.STATVAR.Sout = ground.PARA.albedo .*  forcing.TEMP.Sin;
            ground.STATVAR.Qh = Q_h(ground, forcing);
            ground.STATVAR.Qe_pot = Q_eq_potET(ground, forcing);

            ground = calculateET(ground);
            
            ground.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
        end
        
        function ground = calculateET(ground)
            if ground.STATVAR.Qe_pot > 0
                fraction_T = getET_fraction(ground);
                fraction_E = getET_fraction(ground);
                
                                
                depth_weighting_E = exp(-1./ground.PARA.evaporationDepth .* cumsum(ground.STATVAR.layerThick));  %exponential decrease with depth
                depth_weighting_E(depth_weighting_E<0.05) = 0;
                depth_weighting_E = depth_weighting_E .* ground.STATVAR.layerThick ./ sum(depth_weighting_E .* ground.STATVAR.layerThick,1); %normalize
                               
                depth_weighting_T = exp(-1./ground.PARA.rootDepth .* cumsum(ground.STATVAR.layerThick));
                depth_weighting_T(depth_weighting_T<0.05) = 0;
                depth_weighting_T = depth_weighting_T .* ground.STATVAR.layerThick ./ sum(depth_weighting_T .* ground.STATVAR.layerThick,1);

                fraction_ET = fraction_T .* depth_weighting_T .* ground.PARA.ratioET + fraction_E .* depth_weighting_E .* (1-ground.PARA.ratioET);
                
                ground.STATVAR.Qe = sum(fraction_ET, 1) .* ground.STATVAR.Qe_pot;
                
                fraction_ET = fraction_ET./max(1e-12, sum(fraction_ET, 1));
                
                ground.TEMP.d_water_ET = -ground.STATVAR.Qe ./ ground.CONST.L_v .* fraction_ET;    %in m water per sec
                
            else  %condensation
                ground.STATVAR.Qe = ground.STATVAR.Qe_pot;
                
                ground.TEMP.d_water_ET = ground.STATVAR.T.*0;                
                ground.TEMP.d_water_ET(1,1) = -ground.STATVAR.Qe ./ ground.CONST.L_v; %in m water per sec, put everything in uppermost grid cell
            end
        end
        
        function fraction = getET_fraction(ground)
            saturation = ground.STATVAR.water ./ (ground.STATVAR.layerThick - ground.STATVAR.mineral - ground.STATVAR.organic);
            fraction=double(ground.STATVAR.T>0).*(double(saturation >= ground.STATVAR.field_capacity) + double(saturation < ground.STATVAR.field_capacity).*0.25.*(1-cos(pi().*saturation./ground.STATVAR.field_capacity)).^2);
        end
    
    end
    
end
