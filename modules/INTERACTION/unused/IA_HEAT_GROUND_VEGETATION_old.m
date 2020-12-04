classdef IA_HEAT_GROUND_VEGETATION < matlab.mixin.Copyable
    
    %> Simone: Interaction class between Vegetation and Ground
    %> Soil heat flux (W/m^2) up and down
    
    properties
        PREVIOUS
        NEXT
        IA_PARENT_GROUND
        IA_CHILD_SNOW
        FRACTIONAL_SNOW_COVER
        STATUS
    end
    
    methods
        
        function  get_boundary_condition_m(ia_heat_ground_vegetation)
            
            ground = ia_heat_ground_vegetation.NEXT;
            vegetation = ia_heat_ground_vegetation.PREVIOUS;
            forcing = vegetation.ForcingV;
            
            if ground.IA_CHILD.STATUS == 0 && forcing.TEMP.snowfall > 0  %zero SWE and snowfall occuring, snow CHILD must be initialized
                ground.IA_CHILD.STATUS = 1;
                ia_heat_ground_snow_vegetation = ground.IA_CHILD;
                snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
                ia_heat_ground_snow_vegetation.FRACTIONAL_SNOW_COVER = 0;
                snow = get_boundary_condition_u_create_CHILD(snow, forcing);
                
                ground = get_boundary_condition_u(ground, forcing); %call the native function for the ground class
                
                vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
                %                 vegetation.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* vegetation.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
                
                % wird in get_derivatives_prognostic(ground) dann ground übergeben, hier in fluxes geschrieben um der ordnung zu folgen und negativ
                ground.TEMP.F_ub_water = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(1) + ground.TEMP.F_ub_water;
                ground.TEMP.F_m_water(2:3) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(2)/2;
                ground.TEMP.F_m_water(4:9) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(3)/6;
                ground.TEMP.F_m_water(10:end) = 0;
                
                vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                
                % Simone: Set ground.STATVAR.vegetation.soilvar.t_soisno to actual CG soil temperature
                vegetation.STATVAR.vegetation.soilvar.t_soisno(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.t_soisno(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.t_soisno(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                
            elseif ground.IA_CHILD.STATUS == 2 %non-zero SWE, but snow is still a child
                ia_heat_ground_snow_vegetation = ground.IA_CHILD;
                snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
                
                ground.IA_CHILD.FRACTIONAL_SNOW_COVER = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.STATVAR.ice ./ (ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PARA.swe_per_cell./2);
                snow.IA_PARENT.FRACTIONAL_SNOW_COVER = ground.IA_CHILD.FRACTIONAL_SNOW_COVER;
                fraction_snow = ground.IA_CHILD.FRACTIONAL_SNOW_COVER; %reassign to new variable for code clarity
                
                %call the surface energy balance for both classes
                snow.STATVAR.Lstar = ground.STATVAR.Lstar;  %assign L_star from parent
                
                %partition rainfall
                total_rain = forcing.TEMP.rainfall; % ./ 1000 ./ (24*3600);
                
                forcing.TEMP.rainfall = fraction_snow .* total_rain;
                snow = get_boundary_condition_u_CHILD(snow, forcing); %call the native function for the snow class
                
                forcing.TEMP.rainfall = (1 - fraction_snow) .* total_rain;
                ground = get_boundary_condition_u(ground, forcing); %#ok<*MCSCM> %call the native function for the ground class
                
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
                
%                 % Simone!!
%                 snow.TEMP.F_lb = -ground.TEMP.F_ub;
                
                vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
                %                 vegetation.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* vegetation.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
                
                % wird in get_derivatives_prognostic(ground) dann ground übergeben, hier in fluxes geschrieben um der ordnung zu folgen und negativ
                ground.TEMP.F_ub_water = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(1) + ground.TEMP.F_ub_water;
                ground.TEMP.F_m_water(2:3) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(2)/2;
                ground.TEMP.F_m_water(4:9) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(3)/6;
                ground.TEMP.F_m_water(10:end) = 0;
                
                vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                
                % Simone: Set ground.STATVAR.vegetation.soilvar.t_soisno to actual CG soil temperature
                vegetation.STATVAR.vegetation.soilvar.t_soisno(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.t_soisno(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.t_soisno(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                
                % does IA need output? Simone!!
                %                 ia_heat_ground_vegetation.IA_CHILD_SNOW = snow;
                
            else
                %                 Simone!!
                ground = get_boundary_condition_u(ground, forcing); %call the native function for the ground class
                
                vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
                vegetation.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* vegetation.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
                
                % wird in get_derivatives_prognostic(ground) dann ground übergeben, hier in fluxes geschrieben um der ordnung zu folgen und negativ
                ground.TEMP.F_ub_water = -vegetation.STATVAR.vegetation.mlcanopyinst.etsoi-vegetation.STATVAR.vegetation.soilvar.transp_per_layer(1) + ground.TEMP.F_ub_water;
                ground.TEMP.F_m_water(2:3) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(2)/2;
                ground.TEMP.F_m_water(4:9) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(3)/6;
                ground.TEMP.F_m_water(10:end) = 0;
                
                vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                
                % Simone: Set ground.STATVAR.vegetation.soilvar.t_soisno to actual CG soil temperature
                vegetation.STATVAR.vegetation.soilvar.t_soisno(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.t_soisno(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.t_soisno(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                
            end
        end
        
        %
        %         function ia_heat_ground_vegetation = get_derivative_energy(ia_heat_ground_vegetation)
        %             if ia_heat_ground_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
        %                 snow = ia_heat_ground_vegetation.IA_CHILD_SNOW;
        %                 ground = ia_heat_ground_vegetation.IA_PARENT_GROUND;
        %
        %                 snow = get_derivatives_prognostic_CHILD(snow);
        %                 ground = get_derivative_energy(ground); %call normal function for ground
        %             else
        %                 ia_heat_ground_vegetation.IA_PARENT_GROUND = get_derivative_energy(ia_heat_ground_vegetation.IA_PARENT_GROUND); %call normal function
        %             end
        %         end
        %
        %
        %         function timestep = get_timestep(ia_heat_ground_vegetation)
        %             if ia_heat_ground_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
        %                 timestep_snow = get_timestep_CHILD(ia_heat_ground_vegetation.IA_CHILD_SNOW);
        %                 timestep_ground =  get_timestep(ia_heat_ground_vegetation.IA_PARENT_GROUND);
        %                 timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
        %             else
        %                 timestep =  get_timestep(ia_heat_ground_vegetation.IA_PARENT_GROUND);
        %             end
        %
        %         end
        %
        %         function ia_heat_ground_vegetation = advance_prognostic(ia_heat_ground_vegetation, timestep)
        %             if ia_heat_ground_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
        %                 snow = ia_heat_ground_vegetation.IA_CHILD_SNOW;
        %                 ground = ia_heat_ground_vegetation.IA_PARENT_GROUND;
        %                 snow = advance_prognostic_CHILD(snow, timestep);
        %                 ground =  advance_prognostic(ground, timestep);
        %
        %             else
        %                 ia_heat_ground_vegetation.IA_PARENT_GROUND = advance_prognostic(ia_heat_ground_vegetation.IA_PARENT_GROUND, timestep);
        %             end
        %
        %         end
        %
        %
        %         function ia_heat_ground_vegetation = compute_diagnostic(ia_heat_ground_vegetation)
        %             if ia_heat_ground_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
        %                 snow = ia_heat_ground_vegetation.IA_CHILD_SNOW;
        %                 ground = ia_heat_ground_vegetation.IA_PARENT_GROUND;
        %
        %                 snow = compute_diagnostic_CHILD(snow, forcing);
        %                 ground = compute_diagostic(ground, forcing);
        %
        %                 if snow.STATVAR.waterIce == 0
        %                     ia_heat_ground_vegetation.STATUS = 0;
        %                     snow = initialize_zero_snow(snow, ground); %set all variables to zero
        %                 end
        %             end
        %         end
        %
        %         function ia_heat_ground_vegetation = mix_conductivity(ia_heat_ground_vegetation)
        %             if ia_heat_ground_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
        %                 snow = ia_heat_ground_vegetation.IA_CHILD_SNOW;
        %                 ground = ia_heat_ground_vegetation.IA_PARENT_GROUND;
        %
        %                 snow = conductivity(snow);
        %                 mixed_cond = (snow.STATVAR.layerThick + ground.STATVAR.layerThick(1,1)) .* snow.STATVAR.thermCond .* ground.STATVAR.thermCond(1,1) ./ ...
        %                     (snow.STATVAR.layerThick .* ground.STATVAR.thermCond(1,1) + ground.STATVAR.layerThick(1,1) .* snow.STATVAR.thermCond);
        %                 ground.STATVAR.thermCond(1,1) = mixed_cond;
        %             end
        %         end
        %
        %         function ia_heat_ground_vegetation = check_trigger(ia_heat_ground_vegetation)
        %             if ia_heat_ground_vegetation.STATUS == 1 && ia_heat_ground_vegetation.IA_CHILD_SNOW.STATVAR.ice > ia_heat_ground_vegetation.IA_CHILD_SNOW.PARA.swe_per_cell./2
        %
        %                 ia_heat_ground_vegetation.IA_CHILD_SNOW.PREVIOUS = ia_heat_ground_vegetation.IA_PARENT_GROUND.PREVIOUS;
        %                 ia_heat_ground_vegetation.IA_CHILD_SNOW.PREVIOUS.NEXT = ia_heat_ground_vegetation.IA_CHILD_SNOW;
        %
        %                 ia_heat_ground_vegetation.IA_CHILD_SNOW.NEXT = ia_heat_ground_vegetation.IA_PARENT_GROUND;
        %                 ia_heat_ground_vegetation.IA_PARENT_GROUND.PREVIOUS = ia_heat_ground_vegetation.IA_CHILD_SNOW;
        %
        %
        %                 temp_store = ia_heat_ground_vegetation.IA_CHILD_SNOW; %the snow class
        %
        %                 ia_class = get_IA_class(class(temp_store), class(temp_store.NEXT)); %put snow on top
        %                 temp_store.IA_NEXT = ia_class;
        %                 temp_store.IA_NEXT.PREVIOUS = temp_store;
        %                 temp_store.IA_NEXT.NEXT = temp_store.NEXT;
        %                 temp_store.NEXT.IA_PREVIOUS = temp_store.IA_NEXT;
        %
        %                 ia_heat_ground_vegetation = [];  % destroy the class
        %             end
        %         end
        %
        %     end
        
    end
end
