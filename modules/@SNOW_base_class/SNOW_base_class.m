%base class for a SNOW object with a free water freeze curve and basic function for adjusting the grid- upper
%boundary condition is not specified; %superclass that cannot be run alone



classdef SNOW_base_class < matlab.mixin.Copyable 
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

        %mandatory functions 

       function snow = provide_variables(snow)  %initializes the subvariables as empty arrays 
            snow = provide_PARA(snow);
            snow = provide_CONST(snow);
            snow = provide_STATVAR(snow);
        end
        
        function variable = initialize_from_file(snow, variable, section)
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
            %empty here
        end
        
        function snow = initialize_zero_snow(snow, parentGround)
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
            %empty here, define in subclass
        end
        
        function snow = get_boundary_condition_l(snow, forcing) 
            snow = get_heatFlux_lb(snow, forcing);
        end
        
        
        function snow = get_derivatives_prognostic(snow)
            snow = get_derivative_energy(snow);
        end
        
        function timestep = get_timestep(snow)  %could involve check for several state variables
             timestep = snow.PARA.dE_max ./ (max(abs(snow.TEMP.d_energy) ./ snow.STATVAR.layerThick));
        end
        
        function snow = advance_prognostic(snow, timestep) %real timestep derived as minimum of several classes
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* snow.TEMP.d_energy;
        end
        
        function snow = compute_diagnostic_first_cell(snow, forcing);
            %empty here
        end
        
        function snow = compute_diagnostic(snow, forcing)
            snow = get_T_water(snow); 
            snow = modify_grid(snow);
            snow = conductivity(snow);
        end
        
        
        
        %non-mandatory functions that are required in subclasses
        
        function snow = conductivity(snow)
            snow = conductivity_snow_Yen(snow);
        end
        
        
        function snow = get_derivative_energy(snow)
            fluxes = (snow.STATVAR.T(1:end-1) - snow.STATVAR.T(2:end)) .* snow.STATVAR.thermCond(1:end-1) .* snow.STATVAR.thermCond(2:end) ./...
                (snow.STATVAR.thermCond(1:end-1).* snow.STATVAR.layerThick(2:end)./2 +  snow.STATVAR.thermCond(2:end).* snow.STATVAR.layerThick(1:end-1)./2 );
            
            d_energy=snow.STATVAR.energy.*0;
            d_energy(1) = d_energy(1) + snow.TEMP.heatFlux_ub;
            d_energy(end) = d_energy(end) + snow.TEMP.heatFlux_lb;
            
            if ~isempty(fluxes)
                d_energy(1:end-1) = d_energy(1:end-1) - fluxes(1:end);
                d_energy(2:end) = d_energy(2:end) + fluxes(1:end);
            end
            snow.TEMP.d_energy = d_energy;
        end
        
        
        function snow = get_T_water(snow)
       
            E_frozen = - snow.STATVAR.waterIce .* snow.CONST.L_f;
            
            snow.STATVAR.T = double(snow.STATVAR.energy<=E_frozen) .* (snow.STATVAR.energy-E_frozen) ./ (snow.STATVAR.waterIce .* snow.CONST.c_i);
            snow.STATVAR.water = double(snow.STATVAR.energy > E_frozen) .*  snow.STATVAR.waterIce .* (E_frozen - snow.STATVAR.energy) ./E_frozen;
            snow.STATVAR.ice = double(snow.STATVAR.energy > E_frozen) .*  snow.STATVAR.waterIce .* (snow.STATVAR.energy) ./E_frozen + double(snow.STATVAR.energy <= E_frozen) .* snow.STATVAR.waterIce;
            
            %subtract water----------------
            snow.STATVAR.layerThick = min(snow.STATVAR.layerThick, snow.STATVAR.ice ./ snow.STATVAR.target_density); %adjust so that old density is maintained; do not increase layerThick (when water refreezes)
            max_water = snow.PARA.field_capacity .* (snow.STATVAR.layerThick - snow.STATVAR.ice);
            
            water_left = min(snow.STATVAR.water, max_water);
            excess_water= max(0, snow.STATVAR.water - water_left);  %should be routed later! Here, just added to water_reservoir.
            
            snow.STATVAR.water = water_left;
            
            snow.STATVAR.waterIce = snow.STATVAR.ice + snow.STATVAR.water;
            snow.STATVAR.water_reservoir = snow.STATVAR.water_reservoir + sum(excess_water);
        end

        
        function snow = modify_grid(snow)
            
            if sum(double(snow.STATVAR.ice < 0.5.*snow.PARA.swe_per_cell)) > 0 %reduce
                i=1;
                while i<size(snow.STATVAR.layerThick,1)
                    if snow.STATVAR.ice(i) < 0.5.*snow.PARA.swe_per_cell
                        snow.STATVAR.waterIce(i+1) = snow.STATVAR.waterIce(i+1) + snow.STATVAR.waterIce(i);
                        snow.STATVAR.energy (i+1) = snow.STATVAR.energy (i+1) + snow.STATVAR.energy (i);
                        snow.STATVAR.layerThick(i+1) = snow.STATVAR.layerThick(i+1) + snow.STATVAR.layerThick(i);
                        snow.STATVAR.water(i+1) = snow.STATVAR.water(i+1) + snow.STATVAR.water(i);
                        snow.STATVAR.ice(i+1) = snow.STATVAR.ice(i+1) + snow.STATVAR.ice(i);
                        E_frozen = - snow.STATVAR.waterIce(i+1) .* snow.CONST.L_f;
                        snow.STATVAR.T(i+1) = double(snow.STATVAR.energy(i+1)<=E_frozen) .* (snow.STATVAR.energy(i+1)-E_frozen) ./ (snow.STATVAR.waterIce(i+1) .* snow.CONST.c_i);
                        
                        snow.STATVAR.waterIce(i,:) = [];
                        snow.STATVAR.energy(i,:) = [];
                        snow.STATVAR.layerThick(i,:) = [];
                        snow.STATVAR.water(i,:) = [];
                        snow.STATVAR.ice(i,:) = [];
                        snow.STATVAR.T(i,:) = [];
                    else
                        i=i+1;
                    end
                end
                if size(snow.STATVAR.layerThick,1)>1 && snow.STATVAR.ice(end) < 0.5.*snow.PARA.swe_per_cell   %last cell melts
                    snow.STATVAR.waterIce(end-1) = snow.STATVAR.waterIce(end-1) + snow.STATVAR.waterIce(end);
                    snow.STATVAR.energy (end-1) = snow.STATVAR.energy (end-1) + snow.STATVAR.energy (end);
                    snow.STATVAR.layerThick(end-1) = snow.STATVAR.layerThick(end-1) + snow.STATVAR.layerThick(end);
                    snow.STATVAR.water(end-1) = snow.STATVAR.water(end-1) + snow.STATVAR.water(end);
                    snow.STATVAR.ice(end-1) = snow.STATVAR.ice(end-1) + snow.STATVAR.ice(end);
                    E_frozen = - snow.STATVAR.waterIce(end-1) .* snow.CONST.L_f;
                    snow.STATVAR.T(end-1) = double(snow.STATVAR.energy(end-1)<=E_frozen) .* (snow.STATVAR.energy(end-1)-E_frozen) ./ (snow.STATVAR.waterIce(end-1) .* snow.CONST.c_i);
                    
                    snow.STATVAR.waterIce(end,:) = [];
                    snow.STATVAR.energy(end,:) = [];
                    snow.STATVAR.layerThick(end,:) = [];
                    snow.STATVAR.water(end,:) = [];
                    snow.STATVAR.ice(end,:) = [];
                    snow.STATVAR.T(end,:) = [];
                end
            end
            
            if snow.STATVAR.ice(1) > 1.5.*snow.PARA.swe_per_cell  %expand, check only first cell
                split_fraction = snow.STATVAR.ice(1) ./ snow.PARA.swe_per_cell; %e.g. 1.6
                sf1 = (split_fraction-1)./split_fraction;
                sf2 = 1./split_fraction;
                
                snow.STATVAR.waterIce = [sf1.*snow.STATVAR.waterIce(1); sf2.*snow.STATVAR.waterIce(1); snow.STATVAR.waterIce(2:end)];
                snow.STATVAR.energy = [sf1.*snow.STATVAR.energy(1); sf2.*snow.STATVAR.energy(1); snow.STATVAR.energy(2:end)];
                snow.STATVAR.layerThick = [sf1.*snow.STATVAR.layerThick(1); sf2.*snow.STATVAR.layerThick(1); snow.STATVAR.layerThick(2:end)];
                snow.STATVAR.water = [sf1.*snow.STATVAR.water(1); sf2.*snow.STATVAR.water(1); snow.STATVAR.water(2:end)];
                snow.STATVAR.ice = [sf1.*snow.STATVAR.ice(1); sf2.*snow.STATVAR.ice(1); snow.STATVAR.ice(2:end)];
                snow.STATVAR.T = [snow.STATVAR.T(1); snow.STATVAR.T];
                
            end
        end

    end
end