%> GROUND_SEDIMENT_HEAT_FCURVE_SALTDIFF class for ground with sediment properties,
%> heat and salt transfer, and a freezing charackteristic curve based on measurements
%> 
%> different sediment types can be set using the properties for 
%> soilType, Salinity, fraction of Mineral and Organic Content, Porosity
classdef GROUND_Sediment_heat_fcurve_saltDiff < GROUND_Sediment_heat_fcurve
       % MISSING
    % - initialization in dependence on lat/lon/zsb, forcing, flag for initial
    %   conditions, depth in the stratigraphy
    % - triggered sedimentation/erosion
    
    %     properties %from the superclass
    %         CONST %constants
    %         PARA %external service parameters, all other
    %         STATVAR  %energy, water content, etc.
    %         TEMP  %derivatives in prognostic timestep and optimal timestep
    %                 
    %         PREVIOUS
    %         NEXT
    %         IA_PREVIOUS
    %         IA_NEXT
    %     end

    methods
        %mandatory functions for each class
             
        function ground = provide_variables(ground)  
            % initializes the subvariables as empty arrays
            ground = provide_variables@GROUND_Sediment_heat_fcurve(ground);
            ground = provide_PARA(ground);
            ground = provide_CONST(ground);
            ground = provide_STATVAR(ground);
        end
        
        function variable = initialize_from_file(ground, variable, section)
            % read in parameters from file
            variable = initialize_from_file@GROUND_Sediment_heat_fcurve(ground, variable, section);    
        end
        
        function ground = assign_global_variables(ground, forcing)
            % assigns ground heat flow from the forcing to the ground
            ground = assign_global_variables@GROUND_Sediment_heat_fcurve(ground, forcing);
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            %initialize layerThick, upperPos, lowerPos, etc in the superclass
            ground = initialize_STATVAR_from_file@GROUND_Sediment_heat_fcurve(ground, grid, forcing, depths);

            %overwrite thermal properties with the methods that include
            %salt
            ground = finalize_STATVAR(ground); %assign all variables, that must be calculated or assigned otherwise, including energy, water and ice contents, thermal conductivity
        end

        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            ground = get_boundary_condition_u@GROUND_Sediment_heat_fcurve(ground, forcing);

            %add upper boundary condition for salinity
            ground = get_boundary_salt_u(ground, forcing);
        end

        function ground = get_boundary_condition_l(ground, forcing)
            ground = get_boundary_condition_l@GROUND_Sediment_heat_fcurve(ground, forcing);
            ground.TEMP.saltFlux_lb = 0;
        end


        function ground = get_derivatives_prognostic(ground)
            [divT, divsaltConc] = get_derivative_temperature_salt(ground);
            ground.TEMP.divT = divT;
            ground.TEMP.divsaltConc = divsaltConc;
        end

        function timestep = get_timestep(ground)  %could involve check for several state variables in seconds!

            courant_number_temperature = min(1/2 * ground.STATVAR.c_eff./ground.STATVAR.thermCond(1:end-1) .* (ground.STATVAR.layerThick).^2);
            courant_number_salt = min(1/2 * 1./ground.STATVAR.saltDiff(1:end-1).* (ground.STATVAR.layerThick).^2);

            %make timestep smaller for abrupt changes of salt influx
            %if ground.TEMP.saltFlux_ub ~= 0 %i.e. if we are in a submarine phase
                timestep_min_salt = ground.PARA.dsaltConc_max / max(ground.TEMP.divsaltConc);
                timestep = min([courant_number_temperature, courant_number_salt, timestep_min_salt])/(3600*24); %convert estimate from seconds to days;
            %else
            %   timestep = min(courant_number_temperature, courant_number_salt)/(3600*24); %convert estimate from seconds to days;
            %end

           timestep = max(timestep, 1);%1/(24*60)); %no timesteps below 1/(24*60) day = 1 minute!
     
        end

        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [days] here!
            timestep = timestep*(3600*24); %convert timestep from days to seconds
            %put advancing in time here
            ground.STATVAR.T = ground.STATVAR.T + timestep .* ground.TEMP.divT;
            ground.STATVAR.saltConc = ground.STATVAR.saltConc + timestep .* ground.TEMP.divsaltConc;

        end

        function ground = compute_diagnostic_first_cell(ground, forcing)
            %put stuff that happens only if this cell is the first cell
            %here
            %marine sedimentation could be happening here

            %What needs to happen here?

            %ground = L_star(ground, forcing);
        end

        function ground = compute_diagnostic(ground, forcing)
            %put stuff that happens in every cell here

            %update conductivity, heat capacity and liquid water content
            ground = getThermalProps_wSalt(ground);

            %isostatic movement could happen here

            %What needs to happen here?
        end


    end

end
