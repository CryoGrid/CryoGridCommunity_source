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
            
            % IA class between Vegetation and Ground. 
            %
            %               VEGETATION
            %                   -
            %                     -  
            %                       -
            %                         -                               SNOW -  -  -  - IA_HEAT_GROUND_SNOW_VEGETATION 
            %                           -                           -                  -
            %                             -                       -                -
            %             IA_HEAT_GROUND_VEGETATION       IA_CHILD             -
            %                       -                         -            -        
            %                         -                    -           -
            %                           -                -         - 
            %                             -           -        -
            %                                 GROUND -  -  - 
            
            
            % Checks for Snow status. If Status = -1 Snow exists and covers the ground. Then this IA class does not do anything, because flux between Ground and Snow 
            % and Ground/Snow/Vegetation are calculated in the new Ground/Snow IA IA_HEAT_GROUND_SNOW_VEGETATION which is activeated as soon as Snow is -1. 
            
            
            ground = ia_heat_ground_vegetation.NEXT;
            vegetation = ia_heat_ground_vegetation.PREVIOUS;
            forcing = vegetation.ForcingV;
            
            if ground.IA_CHILD.STATUS == 0 && forcing.TEMP.snowfall > 0 %forcing.TEMP.snowfall > 10 %1.000e-12; 1.000e-3 %instead of 0!!!!!!!  %zero SWE and snowfall occuring, snow CHILD must be initialized
                ground.IA_CHILD.STATUS = 1;
                ia_heat_ground_snow_vegetation = ground.IA_CHILD;
                snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
                ia_heat_ground_snow_vegetation.FRACTIONAL_SNOW_COVER = 0;
 
                forcing.TEMP.rainfall = forcing.TEMP.rainfall + ground.STATVAR.surface_runoff;
                
                snow = get_boundary_condition_u_create_CHILD(snow, forcing);
                
                ground = get_boundary_condition_u(ground, forcing); %call the native function for the ground class
                
                % fluxes from ground and snow overwritten by fluxes from vegetation module
                ground.STATVAR.Qh = vegetation.STATVAR.vegetation.mlcanopyinst.shsoi;
                ground.STATVAR.Qe = vegetation.STATVAR.vegetation.mlcanopyinst.lhsoi;
                ground.TEMP.F_ub = vegetation.STATVAR.vegetation.mlcanopyinst.gsoi;    % forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
%                 ground.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
                
                
                vegetation.STATVAR.vegetation.soilvar.t_top_surfacecell = ground.STATVAR.T(1) + 273.15;
                vegetation.STATVAR.vegetation.soilvar.dz_topsurfacecell = ground.STATVAR.layerThick(1);
                vegetation.STATVAR.vegetation.soilvar.thk_topsurfacecell = ground.STATVAR.thermCond(1);

                vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
                %                 vegetation.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* vegetation.STATVAR.vegetation.mlcanopyinst.etsoi)./1000; %evaporation
                
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
                
                vegetation.STATVAR.vegetation.soilvar.thk(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.thk(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.thk(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                
            elseif ground.IA_CHILD.STATUS == 2  %non-zero SWE, but snow is still a child
                ia_heat_ground_snow_vegetation = ground.IA_CHILD;
                snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
                
                ground.IA_CHILD.FRACTIONAL_SNOW_COVER = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.STATVAR.ice ./ (ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PARA.swe_per_cell./2);
                snow.IA_PARENT.FRACTIONAL_SNOW_COVER = ground.IA_CHILD.FRACTIONAL_SNOW_COVER;
                fraction_snow = ground.IA_CHILD.FRACTIONAL_SNOW_COVER; %reassign to new variable for code clarity
                
                %call the surface energy balance for both classes
                snow.STATVAR.Lstar = ground.STATVAR.Lstar;  %assign L_star from parent
                
                 %partition rainfall %%Simone: add surface runoff
                total_rain = forcing.TEMP.rainfall + ground.STATVAR.surface_runoff + snow.STATVAR.water_reservoir;
 
                forcing.TEMP.rainfall = fraction_snow .* total_rain;
                snow = get_boundary_condition_u_CHILD(snow, forcing); %call the native function for the snow class
               
                forcing.TEMP.rainfall = (1 - fraction_snow) .* total_rain + snow.STATVAR.water_reservoir;

                % fluxes from ground and snow overwritten by fluxes from vegetation module
                snow.STATVAR.Qh = vegetation.STATVAR.vegetation.mlcanopyinst.shsoi;
                snow.STATVAR.Qe = vegetation.STATVAR.vegetation.mlcanopyinst.lhsoi;
                
%                 snow.TEMP.F_ub = vegetation.STATVAR.vegetation.mlcanopyinst.gsoi;     
                snow.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Sout - snow.STATVAR.Qh - snow.STATVAR.Qe;
                                
                forcing.TEMP.rainfall = (1 - fraction_snow) .* total_rain;
                ground = get_boundary_condition_u(ground, forcing); %#call the native function for the ground class
                
                % fluxes from ground and snow overwritten by fluxes from vegetation module
                ground.STATVAR.Qh = vegetation.STATVAR.vegetation.mlcanopyinst.shsoi;
                ground.STATVAR.Qe = vegetation.STATVAR.vegetation.mlcanopyinst.lhsoi;
                ground.TEMP.F_ub = vegetation.STATVAR.vegetation.mlcanopyinst.gsoi;    % forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
%                 ground.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
                
                forcing.TEMP.rainfall = total_rain; %reassign the total rainfall, in case it is used later
                
                %mix the output and assign to the ground class, which calculates L_star later
                ground.STATVAR.Lout = (1-fraction_snow) .* ground.STATVAR.Lout + fraction_snow.* snow.STATVAR.Lout; %mix the surface heat fluxes from snow and ground
                ground.STATVAR.Sout = (1-fraction_snow) .* ground.STATVAR.Sout + fraction_snow.* snow.STATVAR.Sout;
                ground.STATVAR.Qh = (1-fraction_snow) .* ground.STATVAR.Qh + fraction_snow.* snow.STATVAR.Qh;
                ground.STATVAR.Qe = (1-fraction_snow) .* ground.STATVAR.Qe + fraction_snow.* snow.STATVAR.Qe;
                
                %1. conductive fluxes
                flux = (snow.STATVAR.T - ground.STATVAR.T(1)) .* snow.STATVAR.thermCond .* ground.STATVAR.thermCond(1) ./...
                    (snow.STATVAR.thermCond.* ground.STATVAR.layerThick(1)./2 + ground.STATVAR.thermCond(1).* (snow.STATVAR.layerThick./2 ./ fraction_snow) ); %scale the snow cell thickness with fractional cover
                
                snow.TEMP.F_ub = fraction_snow .* (snow.TEMP.F_ub - flux);  %scale the heat flux to a layer with "real layer thickness". This is energy-conserving.
                
                ground.TEMP.F_ub = (1-fraction_snow) .* ground.TEMP.F_ub + fraction_snow .* flux;
                
                vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
                %                 vegetation.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* vegetation.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
                
                vegetation.STATVAR.vegetation.soilvar.t_top_surfacecell = ground.STATVAR.T(1) + 273.15;
                vegetation.STATVAR.vegetation.soilvar.dz_topsurfacecell = ground.STATVAR.layerThick(1);
                vegetation.STATVAR.vegetation.soilvar.thk_topsurfacecell = ground.STATVAR.thermCond(1);

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
                
                vegetation.STATVAR.vegetation.soilvar.thk(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.thk(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.thk(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                
            elseif  ground.IA_CHILD.STATUS == -1 %snow is a child
                % fluxes are calculated in IA_HEAT
                ia_heat_ground_snow_vegetation = ground.IA_CHILD;
                snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
                vegetation = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PREVIOUS;
                forcing = vegetation.ForcingV;
            
%             if ground.IA_CHILD.STATUS == -1 % Child exists -> Ground is snow covered
                
                %%%%%%%%%%%%%% CHANGED HERE: ! %%%%%%%%%%%%%%%
                forcing.TEMP.rainfall = forcing.TEMP.rainfall + ground.STATVAR.surface_runoff;
                snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000./(24.*3600); % Simone!! ./(24.*3600); %snowfall is in mm/day 

                if snow.TEMP.snowfall > 0
                    if forcing.TEMP.snow_reservoir > 1.000e-4 %1.000e-2
                        snow.TEMP.snowfall = snow.TEMP.snowfall + forcing.TEMP.snow_reservoir;
                        forcing.TEMP.snow_reservoir = 0; 
                    end 
                    if snow.TEMP.snowfall < 1.000e-10
                        forcing.TEMP.snow_reservoir = forcing.TEMP.snow_reservoir + snow.TEMP.snowfall;
                        snow.TEMP.snowfall = 0;
                    else
                        snow.TEMP.snowfall = snow.TEMP.snowfall;
                    end
                else
                    snow.TEMP.snowfall = snow.TEMP.snowfall;
                end
                
                vegetation.ForcingV.TEMP.snow_reservoir = forcing.TEMP.snow_reservoir;
                
                disp('snow');
                disp(snow.TEMP.snowfall);
                
                snow = get_boundary_condition_u(snow, forcing);
                snow.STATVAR.Lout = (1-snow.PARA.epsilon) .* forcing.TEMP.Lin + snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+ 273.15).^4;
                snow.STATVAR.Sout = snow.PARA.albedo .*  forcing.TEMP.Sin;
                                                                                         
                snow.TEMP.rainfall = forcing.TEMP.rainfall ./1000./(24.*3600); % Simone!! ./(24.*3600);
                snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);
                snow.TEMP.rain_energy = snow.TEMP.rainfall .* max(0, forcing.TEMP.Tair) .* snow.CONST.c_w;
                
                % fluxes from ground and snow overwritten by fluxes from vegetation module
%                 snow.STATVAR.Qh = vegetation.STATVAR.vegetation.mlcanopyinst.shsoi;
%                 snow.STATVAR.Qe = vegetation.STATVAR.vegetation.mlcanopyinst.lhsoi;
%                 snow.TEMP.F_ub = vegetation.STATVAR.vegetation.mlcanopyinst.gsoi;    % forcing.TEMP.Sin + forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Sout - snow.STATVAR.Qh - snow.STATVAR.Qe;
                snow.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Sout - snow.STATVAR.Qh - snow.STATVAR.Qe;
                
                vegetation.STATVAR.vegetation.soilvar.t_top_surfacecell = snow.STATVAR.T(1) + 273.15;
                vegetation.STATVAR.vegetation.soilvar.dz_topsurfacecell = snow.STATVAR.layerThick(1)/2;
                vegetation.STATVAR.vegetation.soilvar.thk_topsurfacecell = snow.STATVAR.thermCond(1);

                vegetation.STATVAR.vegetation.flux.albsoib = [snow.TEMP.albedo,snow.TEMP.albedo]; % Direct beam albedo of ground (soil)
                vegetation.STATVAR.vegetation.flux.albsoid = [snow.TEMP.albedo,snow.TEMP.albedo]; % Diffuse albedo of ground (soil)
                
                vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
                %                 vegetation.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* vegetation.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
                
                % if waterr < F_ub_water, transp = 0
                ground.TEMP.F_ub_water = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(1); %+ ground.TEMP.F_ub_water;
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
                
                vegetation.STATVAR.vegetation.soilvar.thk(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.thk(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.thk(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
                
            else
                forcing.TEMP.rainfall = forcing.TEMP.rainfall + ground.STATVAR.surface_runoff;

                ground = get_boundary_condition_u(ground, forcing); %call the native function for the ground class
                
                % fluxes from ground and snow overwritten by fluxes from vegetation module
                ground.STATVAR.Qh = vegetation.STATVAR.vegetation.mlcanopyinst.shsoi;
                ground.STATVAR.Qe = vegetation.STATVAR.vegetation.mlcanopyinst.lhsoi;
                ground.TEMP.F_ub = vegetation.STATVAR.vegetation.mlcanopyinst.gsoi;    % forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
%                 ground.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
                
                vegetation.STATVAR.vegetation.soilvar.t_top_surfacecell = ground.STATVAR.T(1) + 273.15;
                vegetation.STATVAR.vegetation.soilvar.dz_topsurfacecell = ground.STATVAR.layerThick(1);
                vegetation.STATVAR.vegetation.soilvar.thk_topsurfacecell = ground.STATVAR.thermCond(1);
                
                vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
                vegetation.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* vegetation.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;

                % wird in get_derivatives_prognostic(ground) dann ground übergeben, hier in fluxes geschrieben um der ordnung zu folgen und negativ
                % vegetation has 3 ground layers
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
                
                vegetation.STATVAR.vegetation.soilvar.thk(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
                vegetation.STATVAR.vegetation.soilvar.thk(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
                vegetation.STATVAR.vegetation.soilvar.thk(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
            end
        end
    end
end
