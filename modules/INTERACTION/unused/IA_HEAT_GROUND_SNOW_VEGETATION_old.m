classdef IA_HEAT_GROUND_SNOW_VEGETATION < matlab.mixin.Copyable
    
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
        
        
        function ia_heat_ground_snow_vegetation = get_boundary_condition_m(ia_heat_ground_snow_vegetation)
           ia_heat_ground_snow_vegetation.IA_PARENT_GROUND = ia_heat_ground_snow_vegetation.NEXT; 
           ground = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;
           ia_heat_ground_snow_vegetation.IA_CHILD_SNOW = ia_heat_ground_snow_vegetation.PREVIOUS;
           snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW; 
            vegetation = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PREVIOUS;
            forcing = vegetation.ForcingV;
%             ia_heat_ground_snow_vegetation.STATUS = 1; % Simone
            
            ia_heat_ground_snow_vegetation.STATUS = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND.IA_CHILD.STATUS;
            
            if (ia_heat_ground_snow_vegetation.STATUS == 1) || (ia_heat_ground_snow_vegetation.STATUS == 2) || (ia_heat_ground_snow_vegetation.STATUS == 0 && forcing.TEMP.snowfall > 0) 
%            if ia_heat_ground_snow_vegetation.STATUS == 1 || forcing.TEMP.snowfall > 0 %non-zero SWE, but snow is still a child, or snowfall occurring
               snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
               ground = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;
               
               ia_heat_ground_snow_vegetation.STATUS = 1;
               ia_heat_ground_snow_vegetation.FRACTIONAL_SNOW_COVER = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.STATVAR.ice ./ (ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PARA.swe_per_cell./2); %
               fraction_snow = ia_heat_ground_snow_vegetation.FRACTIONAL_SNOW_COVER; %reassign to new variable for code clarity
               
               %call the surface energy balance for both classes
               snow.STATVAR.Lstar = ground.STATVAR.Lstar;  %assign L_star from parent
               
               %partition rainfall
%                total_rain = forcing.TEMP.rainfall;
               total_rain_snow = forcing.TEMP.rainfall;
               
               forcing.TEMP.rainfall = fraction_snow .* total_rain_snow;

                
               snow = get_boundary_condition_u(snow, forcing); %call the native function for the ground class
               
%                snow = get_boundary_condition_u(snow, forcing); %call the native function for the snow class   
               
               forcing.TEMP.rainfall = (1 - fraction_snow) .* total_rain_snow;
               ground = get_boundary_condition_u(ground, forcing); %call the native function for the ground class
               
               forcing.TEMP.rainfall = total_rain_snow; %reassign the total rainfall, in case it is used later
               
               %mix the output and assign to the ground class, which calculates L_star later 
               ground.STATVAR.Lout = (1-fraction_snow) .* ground.STATVAR.Lout + fraction_snow.* snow.STATVAR.Lout; %mix the surface heat fluxes from snow and ground
               ground.STATVAR.Sout = (1-fraction_snow) .* ground.STATVAR.Sout + fraction_snow.* snow.STATVAR.Sout; 
               ground.STATVAR.Qh = (1-fraction_snow) .* ground.STATVAR.Qh + fraction_snow.* snow.STATVAR.Qh; 
               ground.STATVAR.Qe = (1-fraction_snow) .* ground.STATVAR.Qe + fraction_snow.* snow.STATVAR.Qe; 
               
               %assign the correct heat fluxes to both classes
               %1. conductive fluxes
               flux = (snow.STATVAR.T - ground.STATVAR.T(1)) .* snow.STATVAR.thermCond .* ground.STATVAR.thermCond(1) ./...
                   (snow.STATVAR.thermCond.* ground.STATVAR.layerThick(1)./2 + ground.STATVAR.thermCond(1).* (snow.STATVAR.layerThick./2 ./ fraction_snow) ); %scale the snow cell thickness with fractional cover
               
               % Simone
               snow.TEMP.F_ub = fraction_snow .* (snow.TEMP.F_ub - flux);  %scale the heat flux to a layer with "real layer thickness". This is energy-conserving.
               snow.TEMP.F_lb = 0;
               ground.TEMP.F_ub = (1-fraction_snow) .* ground.TEMP.F_ub + fraction_snow .* flux;
               
            elseif ia_heat_ground_snow_vegetation.STATUS == -1
               ia_heat_ground_snow_vegetation.IA_CHILD_SNOW = get_boundary_condition_u(ia_heat_ground_snow_vegetation.IA_CHILD_SNOW, forcing); %call the native function for the ground class
            else
               ia_heat_ground_snow_vegetation.IA_PARENT_GROUND = get_boundary_condition_u(ia_heat_ground_snow_vegetation.IA_PARENT_GROUND, forcing); %call the native function for the ground class
           end
        end
        
              function ia_heat_ground_snow_vegetation = get_derivative_energy(ia_heat_ground_snow_vegetation)
            if ia_heat_ground_snow_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
                ground = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;

                snow = get_derivatives_prognostic_CHILD(snow);
                ground = get_derivative_energy(ground); %call normal function for ground
            else
                ia_heat_ground_snow_vegetation.IA_PARENT_GROUND = get_derivative_energy(ia_heat_ground_snow_vegetation.IA_PARENT_GROUND); %call normal function
            end
        end
        
        
        function timestep = get_timestep(ia_heat_ground_snow_vegetation)
            if ia_heat_ground_snow_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
                timestep_snow = get_timestep_CHILD(ia_heat_ground_snow_vegetation.IA_CHILD_SNOW);
                timestep_ground =  get_timestep(ia_heat_ground_snow_vegetation.IA_PARENT_GROUND);
                timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
            else
                timestep =  get_timestep(ia_heat_ground_snow_vegetation.IA_PARENT_GROUND);
            end
           
        end
        
        function ia_heat_ground_snow_vegetation = advance_prognostic(ia_heat_ground_snow_vegetation, timestep)
            if ia_heat_ground_snow_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
                ground = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;
                snow = advance_prognostic_CHILD(snow, timestep);
                ground =  advance_prognostic(ground, timestep);

            else
                ia_heat_ground_snow_vegetation.IA_PARENT_GROUND = advance_prognostic(ia_heat_ground_snow_vegetation.IA_PARENT_GROUND, timestep);
            end
            
        end
        
        
        function ia_heat_ground_snow_vegetation = compute_diagnostic(ia_heat_ground_snow_vegetation)
            if ia_heat_ground_snow_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
                ground = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;
                
                snow = compute_diagnostic_CHILD(snow, forcing);
                ground = compute_diagostic(ground, forcing);
                
                if snow.STATVAR.waterIce == 0
                    ia_heat_ground_snow_vegetation.STATUS = 0;
                    snow = initialize_zero_snow(snow, ground); %set all variables to zero
                end
            end
        end
        
        function ia_heat_ground_snow_vegetation = mix_conductivity(ia_heat_ground_snow_vegetation)
            if ia_heat_ground_snow_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
                snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
                ground = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;
                               
                snow = conductivity(snow);
                mixed_cond = (snow.STATVAR.layerThick + ground.STATVAR.layerThick(1,1)) .* snow.STATVAR.thermCond .* ground.STATVAR.thermCond(1,1) ./ ...
                    (snow.STATVAR.layerThick .* ground.STATVAR.thermCond(1,1) + ground.STATVAR.layerThick(1,1) .* snow.STATVAR.thermCond);
                ground.STATVAR.thermCond(1,1) = mixed_cond;
            end
        end
        
        function ia_heat_ground_snow_vegetation = check_trigger(ia_heat_ground_snow_vegetation)
           if ia_heat_ground_snow_vegetation.STATUS == 1 && ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.STATVAR.ice > ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PARA.swe_per_cell./2
               
               ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PREVIOUS = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND.PREVIOUS;
               ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PREVIOUS.NEXT = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
               
               ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.NEXT = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;
               ia_heat_ground_snow_vegetation.IA_PARENT_GROUND.PREVIOUS = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;

               
               temp_store = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW; %the snow class
               
               ia_class = get_IA_class(class(temp_store), class(temp_store.NEXT)); %put snow on top
               temp_store.IA_NEXT = ia_class;
               temp_store.IA_NEXT.PREVIOUS = temp_store;
               temp_store.IA_NEXT.NEXT = temp_store.NEXT;
               temp_store.NEXT.IA_PREVIOUS = temp_store.IA_NEXT;

               ia_heat_ground_snow_vegetation = [];  % destroy the class
           end
        end
        
    end
    
        
% % %         function  [ia_heat_ground_snow_vegetation] = get_boundary_condition_m(ia_heat_ground_snow_vegetation)
% % %             
% % % %             stratigraphy2 = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
% % % %             stratigraphy1 = ia_heat_ground_snow_vegetation.NEXT;
% % % %             vegetation = stratigraphy2.PREVIOUS;
% % %             
% % %             % forcing should come out of vegetation -> vegetation would be something like ia_heat_ground_snow_vegetation.NEXT
% % %             %    -> rather build a struct in vegetation that can be used here directly ia_heat_ground_snow_vegetation.NEXT.forcing
% % %             %       (from here get_boundary_condition_creat_CHILD can get its forcing)
% % %             % ground comes out of ia_heat_ground_snow_vegetation
% % %             stratigraphy2 = ia_heat_ground_snow_vegetation.PREVIOUS;
% % %             stratigraphy1 = ia_heat_ground_snow_vegetation.NEXT;
% % %             
% % %             vegetation = stratigraphy2.PREVIOUS;
% % %             forcing = vegetation.ForcingV;
% % %             
% % %             if stratigraphy1.IA_CHILD.STATUS == 0 && forcing.TEMP.snowfall > 0  %zero SWE and snowfall occuring, snow CHILD must be initialized
% % %                 stratigraphy1.IA_CHILD.STATUS = 1;
% % %                 ia_heat_ground_snow_vegetation = stratigraphy1.IA_CHILD;
% % %                 stratigraphy2 = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
% % %                 %                 ia_heat_ground_snow_vegetation.FRACTIONAL_SNOW_COVER = 0;
% % %                 stratigraphy2 = get_boundary_condition_u_create_CHILD(stratigraphy2, forcing);
% % %                 
% % %                 % Simone!!
% % %                 stratigraphy1 = get_boundary_condition_u(stratigraphy1, forcing); %call the native function for the ground class
% % %                 
% % %                 %                 % does IA need output? Simone!!
% % %                 %                  ia_heat_ground_snow_vegetation.IA_CHILD_SNOW = stratigraphy2;
% % %                 %                  ia_heat_ground_snow_vegetation.IA_PARENT_GROUND = stratigraphy1;
% % %                 
% % %                 vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
% % %                 %                 stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
% % %                 
% % %                 % wird in get_derivatives_prognostic(ground) dann ground übergeben, hier in fluxes geschrieben um der ordnung zu folgen und negativ
% % %                 stratigraphy2.TEMP.F_ub_water = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(1) + stratigraphy2.TEMP.F_ub_water;
% % %                 stratigraphy2.TEMP.F_m_water(2:3) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(2)/2;
% % %                 stratigraphy2.TEMP.F_m_water(4:9) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(3)/6;
% % %                 stratigraphy2.TEMP.F_m_water(10:end) = 0;
% % %                 
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(1) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(2) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(3) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(1) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(2) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(3) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %                 
% % %                 % Simone: Set ground.STATVAR.vegetation.soilvar.t_soisno to actual CG soil temperature
% % %                 vegetation.STATVAR.vegetation.soilvar.t_soisno(1) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.t_soisno(2) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.t_soisno(3) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %                 
% % %                 %             elseif stratigraphy1.IA_CHILD.STATUS == 2 %non-zero SWE, but snow is still a child
% % %                 %                 ia_heat_ground_snow_vegetation = stratigraphy1.IA_CHILD;
% % %                 %                 stratigraphy2 = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
% % %                 %
% % %                 %                 stratigraphy1.IA_CHILD.FRACTIONAL_SNOW_COVER = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.STATVAR.ice ./ (ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PARA.swe_per_cell./2);
% % %                 %                 stratigraphy2.IA_PARENT.FRACTIONAL_SNOW_COVER = stratigraphy1.IA_CHILD.FRACTIONAL_SNOW_COVER;
% % %                 %                 fraction_snow = stratigraphy1.IA_CHILD.FRACTIONAL_SNOW_COVER; %reassign to new variable for code clarity
% % %                 %
% % %                 %
% % %             elseif stratigraphy1.IA_CHILD.STATUS == 1 || forcing.TEMP.snowfall > 0 %non-zero SWE, but snow is still a child, or snowfall occurring
% % %                 stratigraphy2 = stratigraphy1.IA_CHILD.IA_CHILD_SNOW;
% % %                 stratigraphy1 = stratigraphy1.IA_CHILD.IA_PARENT_GROUND;
% % %                 
% % %                 stratigraphy1.IA_CHILD.STATUS = 1;
% % %                 stratigraphy1.IA_CHILD.FRACTIONAL_SNOW_COVER = stratigraphy1.IA_CHILD.IA_CHILD_SNOW.STATVAR.ice ./ (stratigraphy1.IA_CHILD.IA_CHILD_SNOW.PARA.swe_per_cell./2); %
% % %                 fraction_snow = stratigraphy1.IA_CHILD.FRACTIONAL_SNOW_COVER; %reassign to new variable for code clarity
% % %                 
% % %                 %call the surface energy balance for both classes
% % %                 stratigraphy2.STATVAR.Lstar = stratigraphy1.STATVAR.Lstar;  %assign L_star from parent
% % %                 
% % %                 %partition rainfall
% % %                 total_rain = forcing.TEMP.rainfall; % ./ 1000 ./ (24*3600);
% % %                 
% % %                 forcing.TEMP.rainfall = fraction_snow .* total_rain;
% % %                 stratigraphy2 = get_boundary_condition_u(stratigraphy2, forcing); %call the native function for the snow class
% % %                 
% % %                 forcing.TEMP.rainfall = (1 - fraction_snow) .* total_rain;
% % %                 
% % %                 % Simone!!
% % %                 stratigraphy1 = get_boundary_condition_u(stratigraphy1, forcing); %#ok<*MCSCM> %call the native function for the ground class
% % %                 
% % %                 forcing.TEMP.rainfall = total_rain; %reassign the total rainfall, in case it is used later
% % %                 
% % %                 %mix the output and assign to the ground class, which calculates L_star later
% % %                 stratigraphy1.STATVAR.Lout = (1-fraction_snow) .* stratigraphy1.STATVAR.Lout + fraction_snow.* stratigraphy2.STATVAR.Lout; %mix the surface heat fluxes from snow and ground
% % %                 stratigraphy1.STATVAR.Sout = (1-fraction_snow) .* stratigraphy1.STATVAR.Sout + fraction_snow.* stratigraphy2.STATVAR.Sout;
% % %                 stratigraphy1.STATVAR.Qh = (1-fraction_snow) .* stratigraphy1.STATVAR.Qh + fraction_snow.* stratigraphy2.STATVAR.Qh;
% % %                 stratigraphy1.STATVAR.Qe = (1-fraction_snow) .* stratigraphy1.STATVAR.Qe + fraction_snow.* stratigraphy2.STATVAR.Qe;
% % %                 
% % %                 %                 if abs(stratigraphy1.STATVAR.Qe)>1000 ||  abs(stratigraphy2.STATVAR.Qe)>1000
% % %                 %                     disp(stratigraphy1.STATVAR.Qe)
% % %                 %                     disp(stratigraphy2.STATVAR.Qe)
% % %                 %                     fraction_snow
% % %                 %                     disp(stratigraphy2.STATVAR.T)
% % %                 %                     erffr
% % %                 %                 end
% % %                 
% % %                 %assign the correct heat fluxes to both classes
% % %                 %1. conductive fluxes
% % %                 flux = (stratigraphy2.STATVAR.T - stratigraphy1.STATVAR.T(1)) .* stratigraphy2.STATVAR.thermCond .* stratigraphy1.STATVAR.thermCond(1) ./...
% % %                     (stratigraphy2.STATVAR.thermCond.* stratigraphy1.STATVAR.layerThick(1)./2 + stratigraphy1.STATVAR.thermCond(1).* (stratigraphy2.STATVAR.layerThick./2 ./ fraction_snow) ); %scale the snow cell thickness with fractional cover
% % %                 
% % %                 stratigraphy2.TEMP.F_ub = fraction_snow .* (stratigraphy2.TEMP.F_ub - flux);  %scale the heat flux to a layer with "real layer thickness". This is energy-conserving.
% % %                 stratigraphy2.TEMP.F_lb = 0;
% % %                 stratigraphy1.TEMP.F_ub = (1-fraction_snow) .* stratigraphy1.TEMP.F_ub + fraction_snow .* flux;
% % %                 
% % %                 % Simone!!
% % %                 stratigraphy2.TEMP.F_lb = -stratigraphy1.TEMP.F_ub;
% % %                 vegetation.STATVAR.vegetation.flux.albsoib = [stratigraphy2.TEMP.albedo, stratigraphy2.TEMP.albedo]; % . [0.8,0.8]; % Direct beam albedo of ground (soil)
% % %                 vegetation.STATVAR.vegetation.flux.albsoid = [stratigraphy2.TEMP.albedo, stratigraphy2.TEMP.albedo]; % [0.8,0.8]; % Diffuse albedo of ground (soil)
% % %                 
% % %                 %                 % does IA need output? Simone!!
% % %                 %                  ia_heat_ground_snow_vegetation.IA_CHILD_SNOW = stratigraphy2;
% % %                 %                  ia_heat_ground_snow_vegetation.IA_PARENT_GROUND = stratigraphy1;
% % %                 %
% % %                 vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
% % %                 %                 stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
% % %                 
% % %                 % wird in get_derivatives_prognostic(ground) dann ground übergeben, hier in fluxes geschrieben um der ordnung zu folgen und negativ
% % %                 stratigraphy1.TEMP.F_ub_water = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(1) + stratigraphy1.TEMP.F_ub_water;
% % %                 stratigraphy1.TEMP.F_m_water(2:3) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(2)/2;
% % %                 stratigraphy1.TEMP.F_m_water(4:9) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(3)/6;
% % %                 stratigraphy1.TEMP.F_m_water(10:end) = 0;
% % %                 
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(1) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(2) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(3) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(1) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(2) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(3) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %                 
% % %                 % Simone: Set ground.STATVAR.vegetation.soilvar.t_soisno to actual CG soil temperature
% % %                 vegetation.STATVAR.vegetation.soilvar.t_soisno(1) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.t_soisno(2) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.t_soisno(3) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %                 
% % %                 
% % %             else
% % %                 % Simone!!
% % %                 stratigraphy1 = get_boundary_condition_u(stratigraphy1, forcing); %call the native function for the ground class
% % %                 
% % %                 vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
% % %                 %                 stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* stratigraphy1.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
% % %                 
% % %                 % wird in get_derivatives_prognostic(ground) dann ground übergeben, hier in fluxes geschrieben um der ordnung zu folgen und negativ
% % %                 stratigraphy1.TEMP.F_ub_water = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(1) + stratigraphy1.TEMP.F_ub_water;
% % %                 stratigraphy1.TEMP.F_m_water(2:3) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(2)/2;
% % %                 stratigraphy1.TEMP.F_m_water(4:9) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(3)/6;
% % %                 stratigraphy1.TEMP.F_m_water(10:end) = 0;
% % %                 
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(1) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(2) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(3) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(1) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(2) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(3) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %                 
% % %                 % Simone: Set ground.STATVAR.vegetation.soilvar.t_soisno to actual CG soil temperature
% % %                 vegetation.STATVAR.vegetation.soilvar.t_soisno(1) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.t_soisno(2) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
% % %                 vegetation.STATVAR.vegetation.soilvar.t_soisno(3) = interp1(stratigraphy1.STATVAR.midPoint,stratigraphy1.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
% % %                 
% % %             end
% % %             
% % %             
% % %         end
% % %         
% % %         function ia_heat_ground_snow_vegetation = get_derivative_energy(ia_heat_ground_snow_vegetation)
% % %             if ia_heat_ground_snow_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
% % %                 snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
% % %                 ground = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;
% % %                 
% % %                 snow = get_derivatives_prognostic_CHILD(snow);
% % %                 ground = get_derivative_energy(ground); %call normal function for ground
% % %             else
% % %                 ia_heat_ground_snow_vegetation.IA_PARENT_GROUND = get_derivative_energy(ia_heat_ground_snow_vegetation.IA_PARENT_GROUND); %call normal function
% % %             end
% % %         end
% % %         
% % %         
% % %         function timestep = get_timestep(ia_heat_ground_snow_vegetation)
% % %             if ia_heat_ground_snow_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
% % %                 timestep_snow = get_timestep_CHILD(ia_heat_ground_snow_vegetation.IA_CHILD_SNOW);
% % %                 timestep_ground =  get_timestep(ia_heat_ground_snow_vegetation.IA_PARENT_GROUND);
% % %                 timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
% % %             else
% % %                 timestep =  get_timestep(ia_heat_ground_snow_vegetation.IA_PARENT_GROUND);
% % %             end
% % %             
% % %         end
% % %         
% % %         function ia_heat_ground_snow_vegetation = advance_prognostic(ia_heat_ground_snow_vegetation, timestep)
% % %             if ia_heat_ground_snow_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
% % %                 snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
% % %                 ground = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;
% % %                 snow = advance_prognostic_CHILD(snow, timestep);
% % %                 ground =  advance_prognostic(ground, timestep);
% % %                 
% % %             else
% % %                 ia_heat_ground_snow_vegetation.IA_PARENT_GROUND = advance_prognostic(ia_heat_ground_snow_vegetation.IA_PARENT_GROUND, timestep);
% % %             end
% % %             
% % %         end
% % %         
% % %         
% % %         function ia_heat_ground_snow_vegetation = compute_diagnostic(ia_heat_ground_snow_vegetation)
% % %             if ia_heat_ground_snow_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
% % %                 snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
% % %                 ground = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;
% % %                 
% % %                 snow = compute_diagnostic_CHILD(snow, forcing);
% % %                 ground = compute_diagostic(ground, forcing);
% % %                 
% % %                 if snow.STATVAR.waterIce == 0
% % %                     ia_heat_ground_snow_vegetation.STATUS = 0;
% % %                     snow = initialize_zero_snow(snow, ground); %set all variables to zero
% % %                 end
% % %             end
% % %         end
% % %         
% % %         function ia_heat_ground_snow_vegetation = mix_conductivity(ia_heat_ground_snow_vegetation)
% % %             if ia_heat_ground_snow_vegetation.STATUS == 1 %non-zero SWE, but snow is still a child
% % %                 snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
% % %                 ground = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;
% % %                 
% % %                 snow = conductivity(snow);
% % %                 mixed_cond = (snow.STATVAR.layerThick + ground.STATVAR.layerThick(1,1)) .* snow.STATVAR.thermCond .* ground.STATVAR.thermCond(1,1) ./ ...
% % %                     (snow.STATVAR.layerThick .* ground.STATVAR.thermCond(1,1) + ground.STATVAR.layerThick(1,1) .* snow.STATVAR.thermCond);
% % %                 ground.STATVAR.thermCond(1,1) = mixed_cond;
% % %             end
% % %         end
% % %         
% % %         function ia_heat_ground_snow_vegetation = check_trigger(ia_heat_ground_snow_vegetation)
% % %             if ia_heat_ground_snow_vegetation.STATUS == 1 && ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.STATVAR.ice > ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PARA.swe_per_cell./2
% % %                 
% % %                 ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PREVIOUS = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND.PREVIOUS;
% % %                 ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PREVIOUS.NEXT = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
% % %                 
% % %                 ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.NEXT = ia_heat_ground_snow_vegetation.IA_PARENT_GROUND;
% % %                 ia_heat_ground_snow_vegetation.IA_PARENT_GROUND.PREVIOUS = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
% % %                 
% % %                 
% % %                 temp_store = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW; %the snow class
% % %                 
% % %                 ia_class = get_IA_class(class(temp_store), class(temp_store.NEXT)); %put snow on top
% % %                 temp_store.IA_NEXT = ia_class;
% % %                 temp_store.IA_NEXT.PREVIOUS = temp_store;
% % %                 temp_store.IA_NEXT.NEXT = temp_store.NEXT;
% % %                 temp_store.NEXT.IA_PREVIOUS = temp_store.IA_NEXT;
% % %                 
% % %                 ia_heat_ground_snow_vegetation = [];  % destroy the class
% % %             end
% % %         end
% % %         
% % %     end

end

