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
        
        % 0 no child
        % 1 child initialized / growing
        % 2 Child existing fractionally or melting
        % -1 child existing
        
        
        function ia_heat_ground_snow_vegetation = get_boundary_condition_m(ia_heat_ground_snow_vegetation)
            ground = ia_heat_ground_snow_vegetation.NEXT;
            snow = ia_heat_ground_snow_vegetation.PREVIOUS;
%             snow = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW;
%             vegetation = ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.PREVIOUS;
%             forcing = vegetation.ForcingV;
%             
% %             if ground.IA_CHILD.STATUS == -1 % Child exists -> Ground is snow covered
%                 
%                 %%%%%%%%%%%%%% CHANGED HERE: Overwritting of snow fluxes with vegetation fluxes! %%%%%%%%%%%%%%%
%                 ia_heat_ground_snow_vegetation.IA_CHILD_SNOW = get_boundary_condition_u(ia_heat_ground_snow_vegetation.IA_CHILD_SNOW, forcing);
%                 snow.STATVAR.Lout = (1-snow.PARA.epsilon) .* forcing.TEMP.Lin + snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+ 273.15).^4;
%                 snow.STATVAR.Sout = snow.PARA.albedo .*  forcing.TEMP.Sin;
%                 % forcing.TEMP.snowfall is in mm h2o/day
%                 snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000./(24.*3600); % Simone!! ./(24.*3600); %snowfall is in mm/day 
%                 snow.TEMP.rainfall = forcing.TEMP.rainfall ./1000./(24.*3600); % Simone!! ./(24.*3600);
%                 snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);
%                 snow.TEMP.rain_energy = snow.TEMP.rainfall .* max(0, forcing.TEMP.Tair) .* snow.CONST.c_w;
%                 
%                 % fluxes from ground and snow overwritten by fluxes from vegetation module
%                 snow.STATVAR.Qh = vegetation.STATVAR.vegetation.mlcanopyinst.shsoi;
%                 snow.STATVAR.Qe = vegetation.STATVAR.vegetation.mlcanopyinst.lhsoi;
% %                 snow.TEMP.F_ub = vegetation.STATVAR.vegetation.mlcanopyinst.gsoi;    % forcing.TEMP.Sin + forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Sout - snow.STATVAR.Qh - snow.STATVAR.Qe;
%                 snow.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Sout - snow.STATVAR.Qh - snow.STATVAR.Qe;
%                 
%                 vegetation.STATVAR.vegetation.soilvar.t_top_surfacecell = snow.STATVAR.T(1) + 273.15;
%                 vegetation.STATVAR.vegetation.soilvar.dz_topsurfacecell = snow.STATVAR.layerThick(1)/2;
%                 vegetation.STATVAR.vegetation.soilvar.thk_topsurfacecell = snow.STATVAR.thermCond(1);
%                 %                 vegetation.STATVAR.vegetation.mlcanopyinst.tg_snow = snow.STATVAR.T(1)+273.15;
%                 %                 vegetation.STATVAR.vegetation.mlcanopyinst.snow = 1;
%                 
%                 vegetation.STATVAR.vegetation.flux.albsoib = [ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.TEMP.albedo,ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.TEMP.albedo]; % Direct beam albedo of ground (soil)
%                 vegetation.STATVAR.vegetation.flux.albsoid = [ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.TEMP.albedo,ia_heat_ground_snow_vegetation.IA_CHILD_SNOW.TEMP.albedo]; % Diffuse albedo of ground (soil)
%                 
%                 vegetation.STATVAR.vegetation.soilvar.transp_per_layer = (0.0181528 .* vegetation.STATVAR.vegetation.soilvar.transp_per_layer)./1000;
%                 %                 vegetation.STATVAR.vegetation.mlcanopyinst.etsoi = (0.0181528 .* vegetation.STATVAR.vegetation.mlcanopyinst.etsoi)./1000;
%                 
%                 % if waterr < F_ub_water, transp = 0
%                 ground.TEMP.F_ub_water = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(1); %+ ground.TEMP.F_ub_water;
%                 ground.TEMP.F_m_water(2:3) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(2)/2;
%                 ground.TEMP.F_m_water(4:9) = -vegetation.STATVAR.vegetation.soilvar.transp_per_layer(3)/6;
%                 ground.TEMP.F_m_water(10:end) = 0;
%                 
%                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
%                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
%                 vegetation.STATVAR.vegetation.soilvar.h2osoi_vol(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.water,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
%                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
%                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
%                 vegetation.STATVAR.vegetation.soilvar.h2osoi_ice(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.ice,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
%                 
%                 % Simone: Set ground.STATVAR.vegetation.soilvar.t_soisno to actual CG soil temperature
%                 vegetation.STATVAR.vegetation.soilvar.t_soisno(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
%                 vegetation.STATVAR.vegetation.soilvar.t_soisno(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
%                 vegetation.STATVAR.vegetation.soilvar.t_soisno(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.T+273.15,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
%                 
%                 vegetation.STATVAR.vegetation.soilvar.thk(1) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(1),'nearest');
%                 vegetation.STATVAR.vegetation.soilvar.thk(2) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(2),'nearest');
%                 vegetation.STATVAR.vegetation.soilvar.thk(3) = interp1(ground.STATVAR.midPoint,ground.STATVAR.thermCond,vegetation.STATVAR.vegetation.soilvar.zi(3),'nearest');
%                 
                % This is what happens in IA_HEAT between snow and ground. Above is what should happen between vegetation and snow.
                flux = (snow.STATVAR.T(end) - ground.STATVAR.T(1)) .* snow.STATVAR.thermCond(end) .* ground.STATVAR.thermCond(1) ./...
                    (snow.STATVAR.thermCond(end).* ground.STATVAR.layerThick(1)./2 + ground.STATVAR.thermCond(1).* snow.STATVAR.layerThick(end)./2 );
                
                snow.TEMP.F_lb = -flux;
                ground.TEMP.F_ub = flux;
%             else
%             end
        end
    end
end