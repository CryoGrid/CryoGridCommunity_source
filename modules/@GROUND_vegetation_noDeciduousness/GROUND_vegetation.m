% inherits from GROUND_base_class and adds SEB as upper boundary
% no interaction with snow is possible here

classdef GROUND_vegetation < matlab.mixin.Copyable

    properties
        CONST %constants
        PARA %external service parameters, all other
        STATVAR  %energy, water content, etc.
        ForcingV %forcing variables for ground and snow upper boundary cond
        TEMP  %derivatives in prognostic timestep and optimal timestep
        PREVIOUS
        NEXT
        IA_PREVIOUS
        IA_NEXT
    end
    
    methods
        
        %mandatory functions for each class
        
        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            % % % ground = provide_variables@GROUND_base_class(ground); %call function of the base class
            ground = provide_PARA(ground); %add additional variables
            ground = provide_CONST(ground);
            ground = provide_STATVAR(ground);
            ground = provide_ForcingV(ground);
        end
        
        function variable = initialize_from_file(ground, variable, section)
            class_variables = fieldnames(variable);
            for j=1:size(class_variables,1)
                for k=1:size(section,1)
                    if strcmp(section{k,1}, class_variables(j,1))
                        variable.(class_variables{j}) = section{k,2};
                    end
                end
            end
        end
        
        function ground = assign_global_variables(ground, forcing)
            % % %             ground = assign_global_variables@GROUND_base_class(ground, forcing); %call function of the base class
            ground.PARA.heatFlux_lb = forcing.PARA.heatFlux_lb;
            ground.PARA.airT_height = forcing.PARA.airT_height;
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            % % %             ground = initialize_STATVAR_from_file@GROUND_base_class(ground, grid, forcing, depths);
            variables = fieldnames(ground.STATVAR);
            range = (grid.MIDPOINTS > depths(1,1) & grid.MIDPOINTS <= depths(1,2));
            ground.STATVAR.layerThick = grid.LAYERTHICK(range,1);
            ground.STATVAR.upperPos = forcing.PARA.altitude - depths(1,1);
            ground.STATVAR.lowerPos = forcing.PARA.altitude - depths(1,2);
            
            for j=1:size(variables,1)
                for i=1:size(grid.variable_names,2)
                    if strcmp(variables{j,1}, grid.variable_names{1,i})
                        ground.STATVAR.(variables{j,1}) = grid.variable_gridded(range,i);
                    end
                end
            end
            
            % % %             ground = finalize_STATVAR(ground); %assign all variables, that must be calculated or assigned otherwise, including energy, water and ice contents, thermal conductivity
            ground = set_up_canopy(ground);
            ground = finalize_STATVAR(ground, forcing); %assign all variables, that must be calculated or assigned otherwise
        end
        
        function [ground] = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            
            ground = surface_energy_forest(ground, forcing);
            
            % write these variables out to plot the vegetation on top of ground and snow
            ground.STATVAR.T = (ground.STATVAR.vegetation.mlcanopyinst.tair-273.15)';
            ground.STATVAR.Qe = ground.STATVAR.vegetation.mlcanopyinst.lhsoi; 
            ground.STATVAR.Qh = ground.STATVAR.vegetation.mlcanopyinst.shsoi; 
            ground.STATVAR.layerThick = (ground.STATVAR.vegetation.mlcanopyinst.zw./ground.STATVAR.vegetation.mlcanopyinst.ztop)';
            ground.STATVAR.waterIce = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
            ground.STATVAR.water = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
            ground.STATVAR.ice = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
            
            % Write forcing struct as the input for ground class under vegetation
            ground.ForcingV.TEMP.Tair = ground.STATVAR.vegetation.mlcanopyinst.tveg(1,2)-273.15;
            ground.ForcingV.TEMP.wind = ground.STATVAR.vegetation.mlcanopyinst.wind(1,2);
            ground.ForcingV.TEMP.wind_top = ground.STATVAR.vegetation.mlcanopyinst.wind(1,12); % wind at top of canopy (used for initial snow density calculation)
            ground.ForcingV.TEMP.Sin = ground.STATVAR.vegetation.flux.swdn(1,1) + ground.STATVAR.vegetation.flux.swdn(1,2); %ground.STATVAR.vegetation.flux.swsoi(1,1) + ground.STATVAR.vegetation.flux.swsoi(1,2); % vegetation.mlcanopyinst.sw_prof(1,2,1); %Canopy layer absorbed radiation
            ground.ForcingV.TEMP.Lin = ground.STATVAR.vegetation.flux.irdn; %ground.STATVAR.vegetation.flux.irsoi(1); %vegetation.flux.ir_source(2,1); %Longwave radiation emitted by bottom leaf layer (W/m2)
            ground.ForcingV.TEMP.p = forcing.TEMP.p; % air pressure at reference height (Pa)
            ground.ForcingV.TEMP.snowfall = ground.STATVAR.vegetation.mlcanopyinst.qflx_prec_grnd_snow .* (24 .*3600); % .* (24*3600); qflx_prec_grnd_snow (mm h2o/s) -> .* (24 .*3600) -> mm h2o/day
            ground.ForcingV.TEMP.rainfall = ground.STATVAR.vegetation.mlcanopyinst.qflx_prec_grnd_rain .* (24 .*3600); % .* (24*3600);
            ground.ForcingV.TEMP.q = forcing.TEMP.q; %specific humidity at refernce height (kg/kg)
            ground.ForcingV.TEMP.t = forcing.TEMP.t; % time_snowfall
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
        %             %ground = get_boundary_condition_l@GROUND_base_class(ground, forcing);
        end
        
        
        function ground = get_derivatives_prognostic(ground)
            %ground = get_derivatives_prognostic@GROUND_base_class(ground);
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
            %timestep = get_timestep@GROUND_base_class(ground);
            timestep = 24*60*60;
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!

        end

        function ground = compute_diagnostic_first_cell(ground, forcing)
            %ground = L_star(ground, forcing);
        end
                
        function ground = compute_diagnostic(ground, forcing)

        end
        
        %non-mandatory fucntions
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
        
        function ground = surface_energy_forest(ground, forcing)
            if forcing.TEMP.t >= ground.STATVAR.execution_t || forcing.TEMP.t == forcing.PARA.start_time
                datestr(forcing.TEMP.t)
                
                vegetation = ground.STATVAR.vegetation;
                
                [vegetation] = set_up_forcing(vegetation, forcing);

                [vegetation] = canopy_fluxes_multilayer(vegetation);
                
                ground.STATVAR.vegetation = vegetation;
                
                ground.STATVAR.execution_t = ground.STATVAR.execution_t + 1/24;
       
            end
        end
    end
    
    
    
end
