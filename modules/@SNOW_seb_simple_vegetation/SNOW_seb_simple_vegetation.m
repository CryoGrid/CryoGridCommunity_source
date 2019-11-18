% inhertits from SNOW_base_class and adds SEB as upper boundary
% designed to function as CHILD of a GROUND class that is compatible with
% SNOW classes; compatible with interaction classes IA_SNOW_GROUND and IA_SNOW_GROUND_fcSimple_salt

classdef SNOW_seb_simple_vegetation < SNOW_base_class
    properties
        CONST %constants
        PARA %external service parameters, all other
        STATVAR  %energy, water content, etc.
        TEMP  %derivatives in prognostic timestep and optimal timestep
        NEXT
        PREVIOUS
        IA_NEXT
        IA_PREVIOUS
    end
     
    methods

        %mandatory functions for each class

       function snow = provide_variables(snow)  %initializes the subvariables as empty arrays 
% % %             snow = provide_variables@SNOW_base_class(snow); 
           snow = provide_PARA(snow);
           snow = provide_CONST(snow);
           snow = provide_STATVAR(snow);
        end
        
        function variable = initialize_from_file(snow, variable, section)
            % %             variable = initialize_from_file@SNOW_base_class(snow, variable, section);
                class_variables = fieldnames(variable);
                for j=1:size(class_variables,1)
                    for k=1:size(section,1)
                        if strcmp(section{k,1}, class_variables(j,1))
                            variable.(class_variables{j}) = section{k,2};
                        end
                    end
                end
        end
        
        function snow = assign_global_variables(snow, forcing)
            snow.PARA.airT_height = forcing.PARA.airT_height;
        end
        
        function snow = initialize_zero_snow(snow, parentGround)
% % %             snow = initialize_zero_snow@SNOW_base_class(snow, parentGround);
                snow = provide_STATVAR(snow); %set all variables to empty arrays
                snow.STATVAR.T = 0;
                snow.STATVAR.energy = 0;
                snow.STATVAR.waterIce = 0;
                snow.STATVAR.layerThick = 0;
                snow.STATVAR.ice=0;
                snow.STATVAR.water=0;
                snow.STATVAR.upperPos = parentGround.STATVAR.upperPos;
                snow.STATVAR.lowerPos = parentGround.STATVAR.upperPos;
                snow.STATVAR.water_reservoir = 0;
        end
        
        
        function snow = get_boundary_condition_u(snow, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            
%             disp('Now I get the boundary cond. of snow')

            snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000 ./(24.*3600); %snowfall is in mm/day
            snow.TEMP.rainfall = forcing.TEMP.rainfall ./1000 ./(24.*3600);
            snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);
            snow.TEMP.rain_energy = snow.TEMP.rainfall .* max(0, forcing.TEMP.Tair) .* snow.CONST.c_w;

            % multilayer canopy but with new albedo
        end
        
        function snow = get_boundary_condition_l(snow, forcing)
            % % %             snow = get_boundary_condition_l@SNOW_base_class(snow, forcing);
            snow = get_heatFlux_lb(snow, forcing);
        end
        
        
        function snow = get_derivatives_prognostic(snow)
% % %             snow = get_derivatives_prognostic@SNOW_base_class(snow);
            snow = get_derivative_energy(snow);
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.rain_energy;  %add this here, since it can melt snow and does not change layerThick - must be taken into account for timestep calculation
        end
%         
        function timestep = get_timestep(snow)  %could involve check for several state variables
             timestep1 = get_timestep@SNOW_base_class(snow);
             timestep2 = min((-snow.STATVAR.energy ./ snow.TEMP.d_energy) .*double(snow.TEMP.d_energy>0) + double(snow.TEMP.d_energy<=0).*1e5); %when snow is melting, do not melt more than there is in a grid cell
             timestep = min(timestep1, timestep2);
        end
        
        function snow = advance_prognostic(snow, timestep) %real timestep derived as minimum of several classes
            snow = advance_prognostic@SNOW_base_class(snow, timestep);
            snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + timestep .* (snow.TEMP.snow_energy + snow.TEMP.rain_energy);  %CHANGED!!
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + timestep .* (snow.TEMP.snowfall + snow.TEMP.rainfall);
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + timestep .* snow.TEMP.snowfall ./ (snow.PARA.density ./1000);
            snow.STATVAR.target_density = snow.STATVAR.ice ./ snow.STATVAR.layerThick;
            snow.STATVAR.target_density(1) = (snow.STATVAR.ice(1) + timestep .* snow.TEMP.snowfall) ./ snow.STATVAR.layerThick(1);
        end
        
        function snow = compute_diagnostic_first_cell(snow, forcing)
            snow = L_star(snow, forcing);
        end
        
        function snow = compute_diagnostic(snow, forcing)
            snow = compute_diagnostic@SNOW_base_class(snow, forcing);
            snow = check_trigger(snow);
        end
        
        
        %non_mandatory functions
        
        function snow = surface_energy_balance(snow, forcing)
            snow.STATVAR.Lout = (1-snow.PARA.epsilon) .* forcing.TEMP.Lin + snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+ 273.15).^4;
            snow.STATVAR.Sout = snow.PARA.albedo .*  forcing.TEMP.Sin;
            snow.STATVAR.Qh = Q_h(snow, forcing);
            snow.STATVAR.Qe = Q_eq_potET(snow, forcing);
            
            snow.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Sout - snow.STATVAR.Qh - snow.STATVAR.Qe;
        end
        
        
        function snow = check_trigger(snow)
            if size(snow.STATVAR.energy,1) ==1 && snow.STATVAR.ice < snow.PARA.swe_per_cell./2
                ground = snow.PREVIOUS;
                
                ground.PREVIOUS = []; %reassign ground
                ground.PREVIOUS.NEXT = ground;
                ground.IA_PREVIOUS=[];
                
                snow.NEXT =[]; %reassign snow -> cut all dependencies
                snow.PREVIOUS = ground;
                snow.IA_NEXT =[];
                snow.IA_PREVIOUS =[];
                
                ground.IA_CHILD = IA_SNOW_GROUND_VEGETATION();  %reinitialize interaction class
                ground.IA_CHILD.STATUS = 1; %snow initially active
                ground.IA_CHILD.IA_PARENT_GROUND = ground;  %attach snow and ground to interaction class
                ground.IA_CHILD.IA_CHILD_SNOW = snow;
                
                snow = ground; %assign snow pointer to ground to return to regular stratigraphy
            end
        end


    end
end














