classdef IA_SNOW_GROUND_crocus_SW

    properties
        IA_PARENT_GROUND
        IA_CHILD_SNOW
        FRACTIONAL_SNOW_COVER
        STATUS
    end
    
   methods
        function ia_snow_ground = get_boundary_condition_u(ia_snow_ground, forcing)
           if ia_snow_ground.STATUS == 1 || forcing.TEMP.snowfall > 0 %non-zero SWE, but snow is still a child, or snowfall occurring
               snow = ia_snow_ground.IA_CHILD_SNOW;
               ground = ia_snow_ground.IA_PARENT_GROUND;
               
               ia_snow_ground.STATUS = 1;
               ia_snow_ground.FRACTIONAL_SNOW_COVER = ia_snow_ground.IA_CHILD_SNOW.STATVAR.ice ./ (ia_snow_ground.IA_CHILD_SNOW.PARA.swe_per_cell./2); %
               fraction_snow = ia_snow_ground.FRACTIONAL_SNOW_COVER; %reassign to new variable for code clarity
               
               %call the surface energy balance for both classes
               snow.STATVAR.Lstar = ground.STATVAR.Lstar;  %assign L_star from parent
               
               %partition rainfall
               total_rain_snow = ground.STATVAR.vegetation.mlcanopyinst.qflx_prec_grnd_rain + vegetation.mlcanopyinst.qflx_prec_grnd_snow;
               
               forcing.TEMP.rainfall_snow = fraction_snow .* total_rain_snow;
               
               %change albedo here Simone:
                ia_snow_ground.IA_PARENT_GROUND.STATVAR.vegetation.flux.albsoib = [0.8,0.8]; % Direct beam albedo of ground (soil)
                ia_snow_ground.IA_PARENT_GROUND.STATVAR.vegetation.flux.albsoid = [0.8,0.8]; % Diffuse albedo of ground (soil)
                
                ia_snow_ground.STATUS = 1;
                
                ia_snow_ground.IA_PARENT_GROUND = surface_energy_forest(ia_snow_ground.IA_PARENT_GROUND, forcing);
               
%                snow = get_boundary_condition_u(snow, forcing); %call the native function for the snow class   
               
               forcing.TEMP.rainfall_snow = (1 - fraction_snow) .* total_rain_snow;
               ground = get_boundary_condition_u(ground, forcing); %call the native function for the ground class
               
               forcing.TEMP.rainfall_snow = total_rain_snow; %reassign the total rainfall, in case it is used later
               
               %mix the output and assign to the ground class, which calculates L_star later 
               ground.STATVAR.Lout = (1-fraction_snow) .* ground.STATVAR.Lout + fraction_snow.* snow.STATVAR.Lout; %mix the surface heat fluxes from snow and ground
               ground.STATVAR.Sout = (1-fraction_snow) .* ground.Sout + fraction_snow.* snow.STATVAR.Sout; 
               ground.STATVAR.Qh = (1-fraction_snow) .* ground.STATVAR.Qh + fraction_snow.* snow.STATVAR.Qh; 
               ground.STATVAR.Qe = (1-fraction_snow) .* ground.STATVAR.Qe + fraction_snow.* snow.STATVAR.Qe; 
               
               %assign the correct heat fluxes to both classes
               %1. conductive fluxes
               flux = (snow.STATVAR.T - ground.STATVAR.T(1)) .* snow.STATVAR.thermCond .* ground.STATVAR.thermCond(1) ./...
                   (snow.STATVAR.thermCond.* ground.STATVAR.layerThick(1)./2 + ground.STATVAR.thermCond(1).* (snow.STATVAR.layerThick./2 ./ fraction_snow) ); %scale the snow cell thickness with fractional cover
               
               
               snow.TEMP.F_ub = fraction_snow .* (snow.TEMP.F_ub - flux);  %scale the heat flux to a layer with "real layer thickness". This is energy-conserving.
               snow.TEMP.F_lb = 0;
               ground.TEMP.F_ub = (1-fraction_snow) .* ground.TEMP.F_ub + fraction_snow .* flux;
               
           else
               ia_snow_ground.IA_PARENT_GROUND = get_boundary_condition_u(ia_snow_ground.IA_PARENT_GROUND, forcing); %call the native function for the ground class
           end
        end
        
        
        function ia_snow_ground = get_derivative_energy(ia_snow_ground)
            if ia_snow_ground.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_snow_ground.IA_CHILD_SNOW;
                ground = ia_snow_ground.IA_PARENT_GROUND;

                snow = get_derivatives_prognostic_CHILD(snow);
                ground = get_derivative_energy(ground); %call normal function for ground
            else
                ia_snow_ground.IA_PARENT_GROUND = get_derivative_energy(ia_snow_ground.IA_PARENT_GROUND); %call normal function
            end
        end
        
        
        function timestep = get_timestep(ia_snow_ground)
            if ia_snow_ground.STATUS == 1 %non-zero SWE, but snow is still a child
                timestep_snow = get_timestep_CHILD(ia_snow_ground.IA_CHILD_SNOW);
                timestep_ground =  get_timestep(ia_snow_ground.IA_PARENT_GROUND);
                timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
            else
                timestep =  get_timestep(ia_snow_ground.IA_PARENT_GROUND);
            end
           
        end
        
        function ia_snow_ground = advance_prognostic(ia_snow_ground, timestep)
            if ia_snow_ground.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_snow_ground.IA_CHILD_SNOW;
                ground = ia_snow_ground.IA_PARENT_GROUND;
                snow = advance_prognostic_CHILD(snow, timestep);
                ground =  advance_prognostic(ground, timestep);

            else
                ia_snow_ground.IA_PARENT_GROUND = advance_prognostic(ia_snow_ground.IA_PARENT_GROUND, timestep);
            end
            
        end
        
        
        function ia_snow_ground = compute_diagnostic(ia_snow_ground)
            if ia_snow_ground.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_snow_ground.IA_CHILD_SNOW;
                ground = ia_snow_ground.IA_PARENT_GROUND;
                
                snow = compute_diagnostic_CHILD(snow, forcing);
                ground = compute_diagostic(ground, forcing);
                
                if snow.STATVAR.waterIce == 0
                    ia_snow_ground.STATUS = 0;
                    snow = initialize_zero_snow(snow, ground); %set all variables to zero
                end
            end
        end
        
        function ia_snow_ground = mix_conductivity(ia_snow_ground)
            if ia_snow_ground.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_snow_ground.IA_CHILD_SNOW;
                ground = ia_snow_ground.IA_PARENT_GROUND;
                               
                snow = conductivity(snow);
                mixed_cond = (snow.STATVAR.layerThick + ground.STATVAR.layerThick(1,1)) .* snow.STATVAR.thermCond .* ground.STATVAR.thermCond(1,1) ./ ...
                    (snow.STATVAR.layerThick .* ground.STATVAR.thermCond(1,1) + ground.STATVAR.layerThick(1,1) .* snow.STATVAR.thermCond);
                ground.STATVAR.thermCond(1,1) = mixed_cond;
            end
        end
        
        function ia_snow_ground = check_trigger(ia_snow_ground)
           if ia_snow_ground.STATUS == 1 && ia_snow_ground.IA_CHILD_SNOW.STATVAR.ice > ia_snow_ground.IA_CHILD_SNOW.PARA.swe_per_cell./2
               
               ia_snow_ground.IA_CHILD_SNOW.PREVIOUS = ia_snow_ground.IA_PARENT_GROUND.PREVIOUS;
               ia_snow_ground.IA_CHILD_SNOW.PREVIOUS.NEXT = ia_snow_ground.IA_CHILD_SNOW;
               
               ia_snow_ground.IA_CHILD_SNOW.NEXT = ia_snow_ground.IA_PARENT_GROUND;
               ia_snow_ground.IA_PARENT_GROUND.PREVIOUS = ia_snow_ground.IA_CHILD_SNOW;

               
               temp_store = ia_snow_ground.IA_CHILD_SNOW; %the snow class
               
               ia_class = get_IA_class(class(temp_store), class(temp_store.NEXT)); %put snow on top
               temp_store.IA_NEXT = ia_class;
               temp_store.IA_NEXT.PREVIOUS = temp_store;
               temp_store.IA_NEXT.NEXT = temp_store.NEXT;
               temp_store.NEXT.IA_PREVIOUS = temp_store.IA_NEXT;

               ia_snow_ground = [];  % destroy the class
           end
        end
        
    end
    
end


