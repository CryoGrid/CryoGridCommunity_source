%base class for a GROUND object with a free water freeze curve - upper
%boundary condition is not specified; %superclass that cannot be run alone

classdef GROUND_freezeC_bucketW_seb < matlab.mixin.Copyable
    properties
        CONST %constants
        PARA %external service parameters, all other
        STATVAR  %energy, water content, etc.
        TEMP  %derivatives in prognostic timestep and optimal timestep
        LOOKUP      %lookup tables for conductivity and heat capacity
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
            
            ground = initializeExcessIce(ground); 
            ground = initialize_lookup(ground); % initializes lookup tables for liquid water content, thermal conductivity and heat capacity, as well as vector for T_frozen
            ground.STATVAR.water = ground.STATVAR.waterIce; 
            ground = compute_diagnostic_oldCG(ground); %computes initial values of diagnostic variables
            ground.PARA.airT_height = forcing.PARA.airT_height;
        end       
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            ground = surfaceEnergyBalanceInfiltration(ground, forcing);
            % add rain water
            ground.TEMP.F_ub_water = forcing.TEMP.rainfall./1000 ./ (24.*3600); %in m/sec
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
            ground = get_F_lb(ground, forcing);
            ground.TEMP.F_lb_water = 0;
        end
        
        
        function ground = get_derivatives_prognostic(ground)
                ground = get_derivative_energy(ground);
                ground.TEMP.d_T = ground.TEMP.d_energy ./ ground.STATVAR.heatCapacity ./ ground.STATVAR.layerThick; % derivative of temperature in K/sec 
                ground.TEMP.dwc_dt(1) = ground.TEMP.dwc_dt(1) + ground.TEMP.F_ub_water;  %in m/sec
                ground.TEMP.dwc_dt(end) = ground.TEMP.dwc_dt(end) + ground.TEMP.F_lb_water;  %in m/sec
         end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
            timestep = ground.PARA.dE_max ./ (max(abs(ground.TEMP.d_energy) ./ ground.STATVAR.layerThick));
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            ground.STATVAR.T = ground.STATVAR.T + timestep .* ground.TEMP.d_T;
            % multiply dwc_dt with timestep
            ground.TEMP.dwc = ground.TEMP.dwc_dt .* timestep;  %in m
%             ground.TEMP.dwc_dt = ground.TEMP.dwc_dt.*0; % RBZ 291119, to avoid adding waterfluxes from last snowfree timestep during snow season
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            ground = L_star(ground, forcing);
        end
        
        function ground = compute_diagnostic(ground, forcing)
            %melt Xice and calculate water pool for that melting
            %add to water change per cell, change ground.TEMP.dwc
            
            %do water balance with infiltration
            ground = bucketScheme(ground);
            ground = compute_diagnostic_oldCG(ground);  %conductivity, water content and heat capacity from lookup table
            ground = compute_diagnostic_unfrozenZone(ground); %change conductivity and heat capacity for unfrozen zone
            %recalculate lookup tables when refreezing of cells
            if sum(double(ground.STATVAR.waterIce ./ ground.STATVAR.layerThick ~= ground.LOOKUP.liquidWaterContent(:,end) & ground.STATVAR.T<=0))>0
                disp('infiltration - reinitializing LUT - freezing of infiltrated cell(s)');
                ground = initialize_lookup(ground);
            end
        end
    
        
        
        
        %non-mandatory functions -> required here so that they are usable
        %in subclasses
        

        
        function ground = get_derivative_energy(ground)
            fluxes = (ground.STATVAR.T(1:end-1) - ground.STATVAR.T(2:end)) .* ground.STATVAR.thermCond(1:end-1) .* ground.STATVAR.thermCond(2:end) ./...
                (ground.STATVAR.thermCond(1:end-1).* ground.STATVAR.layerThick(2:end)./2 +  ground.STATVAR.thermCond(2:end).* ground.STATVAR.layerThick(1:end-1)./2 );
            
            d_energy=ground.STATVAR.T.*0;

            d_energy(1) = ground.TEMP.F_ub - fluxes(1);
            d_energy(2:end-1) = fluxes(1:end-1) - fluxes(2:end);
            d_energy(end) = ground.TEMP.F_lb + fluxes(end);
            
            ground.TEMP.d_energy = d_energy;
        end
        
        function ground = get_T_water(ground)  % from new CG without freeze curve, is not used, but required so that the mudule initializes correctly

            Lf = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;

            E_frozen = -Lf.*ground.STATVAR.waterIce;
            
            ground.STATVAR.T = double(ground.STATVAR.energy < E_frozen) .* (ground.STATVAR.energy - E_frozen) ./ (c_i.*ground.STATVAR.waterIce + c_m.*ground.STATVAR.mineral + c_o.*ground.STATVAR.organic) + ...
                double(ground.STATVAR.energy >0) .* ground.STATVAR.energy ./ (c_w.*ground.STATVAR.waterIce + c_m.*ground.STATVAR.mineral + c_o.*ground.STATVAR.organic);
            ground.STATVAR.ice = double(ground.STATVAR.energy <= E_frozen) .*ground.STATVAR.waterIce + double(ground.STATVAR.energy > E_frozen & ground.STATVAR.energy < 0) .* ground.STATVAR.energy ./ (-Lf);
            ground.STATVAR.water = double(ground.STATVAR.energy >= 0) .*ground.STATVAR.waterIce + double(ground.STATVAR.energy > - Lf.*ground.STATVAR.waterIce & ground.STATVAR.energy < 0) .* (ground.STATVAR.energy + Lf.*ground.STATVAR.waterIce) ./ Lf;
            
        end
        
        function ground = conductivity(ground)  % same as above
            ground = conductivity_mixing_squares(ground);
        end
        
    end
    
    
end
