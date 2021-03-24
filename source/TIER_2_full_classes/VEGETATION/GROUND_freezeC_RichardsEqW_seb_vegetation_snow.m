%========================================================================
% CryoGrid GROUND class GROUND_freezeC_RichardsEqW_seb_snow
% heat conduction, Richards equation water scheme, freeze curve based on
% freezing=drying assumption, surface energy balance
% S. Westermann, October 2020
%========================================================================

classdef GROUND_freezeC_RichardsEqW_seb_vegetation_snow < GROUND_freezeC_RichardsEqW_seb_vegetation
    properties
        CHILD
        IA_CHILD
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        

       function ground = provide_PARA(ground)  
            ground = provide_PARA@GROUND_freezeC_RichardsEqW_seb_vegetation(ground);
       end
       
       function ground = provide_CONST(ground)  
           ground = provide_CONST@GROUND_freezeC_RichardsEqW_seb_vegetation(ground);
       end
       
       function ground = provide_STATVAR(ground)  
           ground = provide_STATVAR@GROUND_freezeC_RichardsEqW_seb_vegetation(ground);
       end
       
       function ground = finalize_init(ground, tile)
           ground = finalize_init@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, tile);
           ground.CHILD = 0; % no snow
           ground.IA_CHILD = 0;
       end

       
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)

            if ground.CHILD == 0  %CHILD does not exist
                ground = get_boundary_condition_u@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, tile); %call the native function for the ground class
                
                forcing = ground.VEGETATION.ForcingV;
                
                if forcing.TEMP.snowfall > 0  %create CHILD 
                    CURRENT = ground.PREVIOUS;  %go to Top() and get the stored SNOW class
                    while ~strcmp(class(CURRENT), 'Top')
                        CURRENT = CURRENT.PREVIOUS;
                    end
                    ground.CHILD = copy(CURRENT.STORE.SNOW);
                    ground.CHILD.PARENT = ground;
                    ground.CHILD.NEXT = ground; 
                    ground.IA_CHILD = get_IA_class(class(ground.CHILD), class(ground));
                    ground.IA_CHILD.PREVIOUS = ground.CHILD;
                    ground.IA_CHILD.NEXT = ground; %SNOW CHILD created
                    
                    ground.CHILD = get_boundary_condition_u_create_CHILD(ground.CHILD, tile);  %initialize with fresh snowfall
                end
            else %CHILD exists
                %vegetation part
                frac_snow = ground.CHILD.STATVAR.area ./  ground.STATVAR.area(1,1);
                ground.STATVAR.albedo4vegetation = (1-frac_snow) .* ground.PARA.albedo + frac_snow .* ground.CHILD.STATVAR.albedo;
                ground.STATVAR.emissivity4vegetation = (1-frac_snow) .* ground.PARA.epsilon + frac_snow .* ground.CHILD.PARA.epsilon; %required for reflection of Lin
                ground.STATVAR.Lout4vegetation = (1-frac_snow) .* ground.PARA.epsilon .* ground.CONST.sigma .* (ground.STATVAR.T(1)+273.15).^4 + ...
                    frac_snow .* ground.CHILD.PARA.epsilon .* ground.CHILD.CONST.sigma .* (ground.CHILD.STATVAR.T+273.15).^4; %emitted part of Lout, reflected part in the vegetation
                %Qh and Qe are already mixed from previous timestep
                
                ground.VEGETATION = get_boundary_condition_u(ground.VEGETATION, tile);
                
                %transpiration (subtracted from entire cell)
                range = [1:size(ground.VEGETATION.STATVAR.mlcanopyinst.transpiration,2)]';
                transpiration = ground.VEGETATION.STATVAR.mlcanopyinst.transpiration' .* ground.STATVAR.area(range);
                
                transpiration = min(transpiration, ground.STATVAR.water(range) /(3600.*2)); %hard limitation to prevent crashs: only half of available water can evaporate per 1h timestep
                ground.TEMP.d_water_ET(range) = ground.TEMP.d_water_ET(range) - transpiration;
                ground.TEMP.d_water_ET_energy(range) =  ground.TEMP.d_water_ET_energy(range) - transpiration  .* (double(ground.STATVAR.T(range)>=0) .* ground.CONST.c_w .* ground.STATVAR.T(range) + ...
                    double(ground.STATVAR.T(range)<0) .* ground.CONST.c_i .* ground.STATVAR.T(range));
                
                ground.PARA.airT_height = ground.VEGETATION.STATVAR.mlcanopyinst.zs(1,2);                %set reference height to middle of first canopy layer
                ground.CHILD.PARA.airT_height = ground.VEGETATION.STATVAR.mlcanopyinst.zs(1,2); %could make this dependent on snow height
                %vegetation done

                
                total_area = ground.STATVAR.area; %store the total area of the ground
                total_waterIce = ground.STATVAR.waterIce;
                total_water = ground.STATVAR.water;
                total_ice = ground.STATVAR.ice;
                total_mineral = ground.STATVAR.mineral;
                total_organic = ground.STATVAR.organic;
                
                %split up area in snow-covered (CHILD) and snow-free part (PARENT)
                ground.STATVAR.area = ground.STATVAR.area - ground.CHILD.STATVAR.area(1,1); %replace by snow-free area
                reduction = ground.STATVAR.area(1) ./ total_area(1);
                ground.STATVAR.waterIce = ground.STATVAR.waterIce .* reduction;
                ground.STATVAR.water = ground.STATVAR.water .* reduction;
                ground.STATVAR.ice = ground.STATVAR.ice .* reduction;
                ground.STATVAR.mineral = ground.STATVAR.mineral .* reduction;
                ground.STATVAR.organic = ground.STATVAR.organic .* reduction;

                %SNOW
                ground.CHILD.STATVAR.Lstar = ground.STATVAR.Lstar;
                ground.CHILD = get_boundary_condition_u_CHILD(ground.CHILD, tile); %reads correct below-canopy forcing
                
                %GROUND
                forcing = ground.VEGETATION.ForcingV;
                
                ground = calculate_radiation(ground, forcing);
                
                ground = Q_evap_CLM4_5(ground, forcing);
                
                ground.STATVAR.Qh = Q_h(ground, forcing);
            
                %energy
                ground.TEMP.F_ub = (forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe) .* ground.STATVAR.area(1);
                ground.TEMP.d_energy(1) = ground.TEMP.d_energy(1) + ground.TEMP.F_ub;
                
                %water -> evaporation
                ground.TEMP.d_water_ET(1,1) = ground.TEMP.d_water_ET(1,1) -  ground.STATVAR.evap.* ground.STATVAR.area(1,1); %in m3 water per sec, put everything in uppermost grid cell
                ground.TEMP.d_water_ET_energy(1,1) = ground.TEMP.d_water_ET_energy(1,1) -  ground.STATVAR.evap_energy.* ground.STATVAR.area(1,1);
                %mass balance of sublimation not considered (for GROUND)!
                
                %upper BC water (rainfall)
                ground = get_boundary_condition_u_RichardsEq(ground, forcing); %checked that this flux can be taken up!!
                
                get_IA_CHILD_boundary_condition_u(ground.IA_CHILD, tile);
                %call designated mandatory function for CHILD-PARENT interactions in
                %the IA class governing IA between SNOW and GROUND                

                ground.STATVAR.Lout = (ground.STATVAR.area(1,1) .* ground.STATVAR.Lout + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Lout) ./ total_area(1,1); %mix the surface heat fluxes from snow and ground
                ground.STATVAR.Sout = (ground.STATVAR.area(1,1) .* ground.STATVAR.Sout + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Sout) ./ total_area(1,1);
                ground.STATVAR.Qh = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qh + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qh) ./ total_area(1,1);
                ground.STATVAR.Qe = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qe + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qe) ./ total_area(1,1);
                
                %reassign the true totals of ground
                ground.STATVAR.area = total_area;
                ground.STATVAR.waterIce = total_waterIce;
                ground.STATVAR.water = total_water;
                ground.STATVAR.ice = total_ice;
                ground.STATVAR.mineral = total_mineral;
                ground.STATVAR.organic = total_organic;

            end
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
              ground = get_boundary_condition_l@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, tile);
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            if ground.CHILD == 0  
                ground = get_derivatives_prognostic@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, tile); %call normal function
            else
                ground.CHILD = get_derivatives_prognostic_CHILD(ground.CHILD, tile);
                ground = get_derivatives_prognostic@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, tile); 
            end
        end
        
        function timestep = get_timestep(ground, tile) 
            if ground.CHILD == 0
                timestep = get_timestep@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, tile);
            else 
                timestep_snow = get_timestep_CHILD(ground.CHILD, tile);
                timestep_ground =  get_timestep@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, tile);
                timestep = timestep_ground + double(timestep_snow > 0 && timestep_snow < timestep_ground) .* (timestep_snow - timestep_ground);
            end
        end
        
        function ground = advance_prognostic(ground, tile) 
            if ground.CHILD == 0
                ground =  advance_prognostic@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, tile);
            else                
                ground.CHILD = advance_prognostic_CHILD(ground.CHILD, tile);
                ground =  advance_prognostic@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, tile);
            end
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            forcing = tile.FORCING;
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, tile)
            if ground.CHILD == 0
                ground = compute_diagnostic@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, tile);
            else
                            %test with bypass flow

                i=1;
                while ground.CHILD.STATVAR.excessWater > 0
                    water_deficit = ground.STATVAR.layerThick(i,1).* ground.STATVAR.area(i,1) - ground.STATVAR.mineral(i,1) - ground.STATVAR.organic(i,1) - ground.STATVAR.waterIce(i,1);
                    infiltrate = min(water_deficit, ground.CHILD.STATVAR.excessWater);
                    ground.STATVAR.waterIce(i,1) = ground.STATVAR.waterIce(i,1) + infiltrate;
                    ground.CHILD.STATVAR.excessWater = ground.CHILD.STATVAR.excessWater - infiltrate;
                    i=i+1;
                end

                
                ground = compute_diagnostic@GROUND_freezeC_RichardsEqW_seb_vegetation(ground, tile);
                ground.CHILD = compute_diagnostic_CHILD(ground.CHILD, tile);
                
            end
        end
        
        function ground = check_trigger(ground, tile)
            
            
            if ground.CHILD ~= 0
                %delete CHILD
                if ground.CHILD.STATVAR.area ./ ground.STATVAR.area(1,1) < 1e-6 %cutoff to get rid of remaining snow
                   ground.CHILD = 0;
                   ground.IA_CHILD = 0;
                %make SNOW CHILD full class   
                elseif ground.CHILD.STATVAR.area ./ ground.STATVAR.area(1,1) > 1 
                    %transforms dimensions and STAVAR
                    snow_volume = ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.layerThick;
                    ground.CHILD.STATVAR.area = ground.STATVAR.area(1,1);
                    ground.CHILD.STATVAR.layerThick = snow_volume ./ ground.CHILD.STATVAR.area;
                   
                    %make snow a real class
                    %reset vegetation
                    ground.VEGETATION.PARENT_SURFACE = ground.CHILD;
                    ground.CHILD.VEGETATION = ground.VEGETATION;
                    
                    ground.CHILD.PARENT = 0;
                    ground.CHILD.PREVIOUS = ground.PREVIOUS;
                    ground.CHILD.NEXT = ground;
                    ground.PREVIOUS.NEXT = ground.CHILD;
                    ground.PREVIOUS = ground.CHILD;
                    
                    ground.CHILD = 0;
                    ground.IA_PREVIOUS = ground.IA_CHILD; 
                    ground.PREVIOUS.IA_NEXT = ground.IA_CHILD;
                    ground.IA_CHILD = 0;
                    
                    
                end
            end
        end
        
    end
end