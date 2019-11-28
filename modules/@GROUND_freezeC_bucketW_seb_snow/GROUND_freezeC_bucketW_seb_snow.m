%base class for a GROUND object with a free water freeze curve - upper
%boundary condition is not specified; %superclass that cannot be run alone

classdef GROUND_freezeC_bucketW_seb_snow < GROUND_freezeC_bucketW_seb
    properties
        IA_CHILD
    end
    
    
    methods
        
        %mandatory functions for each class
        

       function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            ground = provide_variables@GROUND_freezeC_bucketW_seb(ground);
        end
        
        function variable = initialize_from_file(ground, variable, section)
            variable = initialize_from_file@GROUND_freezeC_bucketW_seb(ground, variable, section);
        end
        
        function ground = assign_global_variables(ground, forcing)
            ground = assign_global_variables@GROUND_freezeC_bucketW_seb(ground, forcing);
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            ground = initialize_STATVAR_from_file@GROUND_freezeC_bucketW_seb(ground, grid, forcing, depths);
        end
        
        function ground = get_boundary_condition_u(ground, forcing)
           
            
            if ground.IA_CHILD.STATUS == 0 && forcing.TEMP.snowfall > 0  %zero SWE and snowfall occuring, snow CHILD must be initialized
               ground.IA_CHILD.STATUS = 1;
               ia_snow_ground = ground.IA_CHILD; 
               snow = ia_snow_ground.IA_CHILD_SNOW;
               ia_snow_ground.FRACTIONAL_SNOW_COVER = 0;
               snow = get_boundary_condition_u_create_CHILD(snow, forcing);
               ground = get_boundary_condition_u@GROUND_freezeC_bucketW_seb(ground, forcing); %call the native function for the ground class
            
            elseif ground.IA_CHILD.STATUS == 2 %non-zero SWE, but snow is still a child
                ia_snow_ground = ground.IA_CHILD; 
                snow = ia_snow_ground.IA_CHILD_SNOW;
                

                ground.IA_CHILD.FRACTIONAL_SNOW_COVER = ia_snow_ground.IA_CHILD_SNOW.STATVAR.ice ./ (ia_snow_ground.IA_CHILD_SNOW.PARA.swe_per_cell./2); 
                snow.IA_PARENT.FRACTIONAL_SNOW_COVER = ground.IA_CHILD.FRACTIONAL_SNOW_COVER;
                fraction_snow = ground.IA_CHILD.FRACTIONAL_SNOW_COVER; %reassign to new variable for code clarity
                
                %call the surface energy balance for both classes
                snow.STATVAR.Lstar = ground.STATVAR.Lstar;  %assign L_star from parent
                
                %partition rainfall
                total_rain = forcing.TEMP.rainfall;
                
                forcing.TEMP.rainfall = fraction_snow .* total_rain;
                snow = get_boundary_condition_u_CHILD(snow, forcing); %call the native function for the snow class
                
                forcing.TEMP.rainfall = (1 - fraction_snow) .* total_rain;
                ground = get_boundary_condition_u@GROUND_freezeC_bucketW_seb(ground, forcing); %call the native function for the ground class
                
                
                forcing.TEMP.rainfall = total_rain; %reassign the total rainfall, in case it is used later
                
                %mix the output and assign to the ground class, which calculates L_star later
                ground.STATVAR.Lout = (1-fraction_snow) .* ground.STATVAR.Lout + fraction_snow.* snow.STATVAR.Lout; %mix the surface heat fluxes from snow and ground
                ground.STATVAR.Sout = (1-fraction_snow) .* ground.STATVAR.Sout + fraction_snow.* snow.STATVAR.Sout;
                ground.STATVAR.Qh = (1-fraction_snow) .* ground.STATVAR.Qh + fraction_snow.* snow.STATVAR.Qh;
                ground.STATVAR.Qe = (1-fraction_snow) .* ground.STATVAR.Qe + fraction_snow.* snow.STATVAR.Qe;
                
%                 if abs(ground.STATVAR.Qe)>1000 ||  abs(snow.STATVAR.Qe)>1000
%                     disp(ground.STATVAR.Qe)
%                     disp(snow.STATVAR.Qe)
%                     fraction_snow
%                     disp(snow.STATVAR.T)
%                     erffr
%                 end
                
                %assign the correct heat fluxes to both classes
                %1. conductive fluxes
                flux = (snow.STATVAR.T - ground.STATVAR.T(1)) .* snow.STATVAR.thermCond .* ground.STATVAR.thermCond(1) ./...
                    (snow.STATVAR.thermCond.* ground.STATVAR.layerThick(1)./2 + ground.STATVAR.thermCond(1).* (snow.STATVAR.layerThick./2 ./ fraction_snow) ); %scale the snow cell thickness with fractional cover
                
                snow.TEMP.F_ub = fraction_snow .* (snow.TEMP.F_ub - flux);  %scale the heat flux to a layer with "real layer thickness". This is energy-conserving.

                ground.TEMP.F_ub = (1-fraction_snow) .* ground.TEMP.F_ub + fraction_snow .* flux;
                
            else
                ground = get_boundary_condition_u@GROUND_freezeC_bucketW_seb(ground, forcing); %call the native function for the ground class
            end
            
            
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
              ground = get_boundary_condition_l@GROUND_freezeC_bucketW_seb(ground, forcing);

        end
        
        
        function ground = get_derivatives_prognostic(ground)
            if ground.IA_CHILD.STATUS <= 1  
                ground = get_derivatives_prognostic@GROUND_freezeC_bucketW_seb(ground); %call normal function
            else
                if ground.IA_CHILD.STATUS == 2 %non-zero SWE, but snow is still a child
                    ia_snow_ground = ground.IA_CHILD;
                    snow = ia_snow_ground.IA_CHILD_SNOW;
                    
                    snow = get_derivatives_prognostic_CHILD(snow);
                    ground = get_derivatives_prognostic@GROUND_freezeC_bucketW_seb(ground); %call normal function for ground
                else
                    ground = get_derivatives_prognostic@GROUND_freezeC_bucketW_seb(ground); %call normal function
                end
            end
        end
        
        function timestep = get_timestep(ground) 
            if ground.IA_CHILD.STATUS <= 1  
                timestep =  get_timestep@GROUND_freezeC_bucketW_seb(ground);
            else  %snow is CHILD
                ia_snow_ground = ground.IA_CHILD;
                snow = ia_snow_ground.IA_CHILD_SNOW;
                
                timestep_snow = get_timestep_CHILD(snow);
                timestep_ground =  get_timestep@GROUND_freezeC_bucketW_seb(ground);
                timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
            end
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            if ground.IA_CHILD.STATUS <= 0
                ground =  advance_prognostic@GROUND_freezeC_bucketW_seb(ground, timestep);
            elseif ground.IA_CHILD.STATUS == 2
                ia_snow_ground = ground.IA_CHILD;
                snow = ia_snow_ground.IA_CHILD_SNOW;
                
                snow = advance_prognostic_CHILD(snow, timestep);
                ground =  advance_prognostic@GROUND_freezeC_bucketW_seb(ground, timestep);
                
            else %ground.IA_CHILD.STATUS == 1
                ia_snow_ground = ground.IA_CHILD;
                snow = ia_snow_ground.IA_CHILD_SNOW;
                
                snow = advance_prognostic_create_CHILD(snow, timestep);
                ground.IA_CHILD.STATUS = 2; %snow exists as normal CHILD from that point on  
            end
        end
    
        
        function ground = compute_diagnostic_first_cell(ground, forcing);
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, forcing)
            if ground.IA_CHILD.STATUS <= 0
                ground = compute_diagnostic@GROUND_freezeC_bucketW_seb(ground, forcing);
            else
                
                ia_snow_ground = ground.IA_CHILD;
                snow = ia_snow_ground.IA_CHILD_SNOW;
                
                ground.IA_CHILD.IA_CHILD_SNOW = compute_diagnostic_CHILD(ground.IA_CHILD.IA_CHILD_SNOW, forcing);
                if snow.STATVAR.water < 0
                    disp(snow.STATVAR.T)
                    disp(snow.STATVAR.energy)
                    disp(snow.STATVAR.waterIce)
                    disp(snow.STATVAR.water)
                    disp(snow.STATVAR.ice)
                    disp(snow.STATVAR.layerThick)
                    disp('next')
                end
                ground = compute_diagnostic@GROUND_freezeC_bucketW_seb(ground, forcing);
            end
        end
        
        
        
        
        %non-mandatory functions -> required here so that they are usable
        %in subclasses
        
        
        
        function ground = get_derivative_energy(ground)
            fluxes = (ground.STATVAR.T(1:end-1) - ground.STATVAR.T(2:end)) .* ground.STATVAR.thermCond(1:end-1) .* ground.STATVAR.thermCond(2:end) ./...
                (ground.STATVAR.thermCond(1:end-1).* ground.STATVAR.layerThick(2:end)./2 +  ground.STATVAR.thermCond(2:end).* ground.STATVAR.layerThick(1:end-1)./2 );
            
            d_energy=ground.STATVAR.T.*0;
            
            d_energy(1) = ground.TEMP.F_ub - fluxes(1);
            d_energy(2:end-1) = fluxes(1:end-1) - fluxes(2:end);
            d_energy(end) = ground.TEMP.F_lb + fluxes(end);
            
            ground.TEMP.d_energy = d_energy;
        end
        
        function ground = get_T_water(ground)  % from new CG without freeze curve, is not used, but required so that the mudule initializes correctly

            Lf = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;

            E_frozen = -Lf.*ground.STATVAR.waterIce;
            
            ground.STATVAR.T = double(ground.STATVAR.energy < E_frozen) .* (ground.STATVAR.energy - E_frozen) ./ (c_i.*ground.STATVAR.waterIce + c_m.*ground.STATVAR.mineral + c_o.*ground.STATVAR.organic) + ...
                double(ground.STATVAR.energy >0) .* ground.STATVAR.energy ./ (c_i.*ground.STATVAR.waterIce + c_m.*ground.STATVAR.mineral + c_o.*ground.STATVAR.organic);
            ground.STATVAR.ice = double(ground.STATVAR.energy <= E_frozen) .*ground.STATVAR.waterIce + double(ground.STATVAR.energy > E_frozen & ground.STATVAR.energy < 0) .* ground.STATVAR.energy ./ (-Lf);
            ground.STATVAR.water = double(ground.STATVAR.energy >= 0) .*ground.STATVAR.waterIce + double(ground.STATVAR.energy > - Lf.*ground.STATVAR.waterIce & ground.STATVAR.energy < 0) .* (ground.STATVAR.energy + Lf.*ground.STATVAR.waterIce) ./ Lf;
            
        end
        
        function ground = conductivity(ground)  % same as above
            ground = conductivity_mixing_squares(ground);
        end
        
    end
    
    
end
