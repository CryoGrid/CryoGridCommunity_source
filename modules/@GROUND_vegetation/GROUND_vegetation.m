% inherits from GROUND_base_class and adds SEB as upper boundary
% no interaction with snow is possible here

classdef GROUND_vegetation < matlab.mixin.Copyable
    
    %%% < GROUND_base_class
    
    properties
        CONST %constants
        PARA %external service parameters, all other
        STATVAR  %energy, water content, etc.
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
            if forcing.TEMP.t >= ground.STATVAR.execution_t
                datestr(forcing.TEMP.t)
                
                vegetation = ground.STATVAR.vegetation;
                
                [vegetation] = set_up_forcing(vegetation, forcing);
                
                [vegetation] = canopy_fluxes_multilayer(vegetation);
                
%                 [vegetation] = figures(vegetation);
                
                ground.STATVAR.vegetation = vegetation;
                ground.STATVAR.execution_t = ground.STATVAR.execution_t + 1/24;
            end
        
            
            %             ground.TEMP.F_ub = ground.STATVAR.vegetation.mlcanopyinst.gsoi;
            %             ground.STATVAR.T(end) = ground.STATVAR.vegetation.mlcanopyinst.tveg(:,2,1)-273.15;
            %             ground.STATVAR.thermCond(end) = 0.025; % Thermal conductivity of air
            %             ground.STATVAR.Lout = ground.STATVAR.vegetation.mlcanopyinst.ircan;
            %             ground.STATVAR.Qh = ground.STATVAR.vegetation.mlcanopyinst.shveg;
            %             ground.STATVAR.Qe = ground.STATVAR.vegetation.mlcanopyinst.lhveg;
        end
        
        %         function ground = get_boundary_condition_l(ground, forcing)
        %             %ground = get_boundary_condition_l@GROUND_base_class(ground, forcing);
        %end
        
        
        function ground = get_derivatives_prognostic(ground)
            %ground = get_derivatives_prognostic@GROUND_base_class(ground);
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
            %timestep = get_timestep@GROUND_base_class(ground);
            timestep = 24*60*60;
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!

            %ground = advance_prognostic@GROUND_base_class(ground, timestep);
            %             ground.STATVAR.current_t = forcing.TEMP.t + timestep;
        end
        
        %         function ground = compute_diagnostic_first_cell(ground, forcing)
        %             %ground = L_star(ground, forcing);
        %         end
        
        %         function ground = compute_diagnostic(ground, forcing)
        %             %ground = compute_diagnostic@GROUND_base_class(ground, forcing);
        %         end
        
        
        %non-mandatory fucntions
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
    end
    
    
    
    
end
