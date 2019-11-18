%base class for a GROUND object with a free water freeze curve - upper
%boundary condition is not specified; %superclass that cannot be run alone

classdef GROUND_base_class < matlab.mixin.Copyable
    
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
            ground = provide_PARA(ground);
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
            ground.PARA.heatFlux_lb = forcing.PARA.heatFlux_lb;
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
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
            ground = finalize_STATVAR(ground); %assign all variables, that must be calculated or assigned otherwise, including energy, water and ice contents, thermal conductivity
        end
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            %empty here
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
            ground = get_F_lb(ground, forcing);
        end
        
        
        function ground = get_derivatives_prognostic(ground)
                ground = get_derivative_energy(ground);
         end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
            timestep = ground.PARA.dE_max ./ (max(abs(ground.TEMP.d_energy) ./ ground.STATVAR.layerThick));
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            ground.STATVAR.energy = ground.STATVAR.energy + timestep .* ground.TEMP.d_energy;
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            %empty here
        end
        
        function ground = compute_diagnostic(ground, forcing)
            ground = get_T_water(ground);
            ground = conductivity(ground);
        end
           
        %non-mandatory functions -> required here so that they are usable
        %in subclasses
        
        function ground = conductivity(ground)
            ground = conductivity_mixing_squares(ground);
        end
        
        function ground = get_derivative_energy(ground)
            fluxes = (ground.STATVAR.T(1:end-1) - ground.STATVAR.T(2:end)) .* ground.STATVAR.thermCond(1:end-1) .* ground.STATVAR.thermCond(2:end) ./...
                (ground.STATVAR.thermCond(1:end-1).* ground.STATVAR.layerThick(2:end)./2 +  ground.STATVAR.thermCond(2:end).* ground.STATVAR.layerThick(1:end-1)./2 );
            
            d_energy=ground.STATVAR.energy.*0;

            d_energy(1) = ground.TEMP.F_ub - fluxes(1);
            d_energy(2:end-1) = fluxes(1:end-1) - fluxes(2:end);
            d_energy(end) = ground.TEMP.F_lb + fluxes(end);
            
            ground.TEMP.d_energy = d_energy;
        end
        
        function ground = get_T_water(ground)

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
        
    end
    
    
end
