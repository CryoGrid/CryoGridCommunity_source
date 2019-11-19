classdef GROUND_subsi < GROUND_base_class
    %GROUND_SUBSIDENCE Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        
        %mandatory functions for each class
        
        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            
            ground = provide_variables@GROUND_base_class(ground); %call function of the base class
            ground = provide_PARA(ground); %add additional variables
            ground = provide_CONST(ground);
            ground = provide_STATVAR(ground);
        end
        
        function variable = initialize_from_file(ground, variable, section, grid)
            variable = initialize_from_file@GROUND_base_class(ground, variable, section);
            %              ground = initializeSoilThermalProperties(ground, grid);
        end
        
        function ground = assign_global_variables(ground, forcing)
            ground = assign_global_variables@GROUND_base_class(ground, forcing);
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            ground = initialize_STATVAR_from_file@GROUND_base_class(ground, grid, forcing, depths);
            ground = initializeSoilThermalProperties(ground, grid); %NC added
        end
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            %assign upper boundary for heat in subclass
            %assign upper boundary for salt in subclass
           ground = getThermalPropertiesInfiltration(ground, forcing);
            
        end
        
        function ground = get_boundary_condition_l(ground)
            ground = get_boundary_condition_l@GROUND_base_class(ground);
            %assign lower boundary for salt in subclass
        end
        
        function ground = get_derivatives_prognostic(ground)
            ground = get_derivative_energy(ground);
            ground = get_derivative_salt(ground);
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables, in this case only check energy derivative
            timestep = get_timestep@GROUND_base_class(ground);
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            ground = advance_prognostic@GROUND_base_class(ground, timestep);
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            %assigned in subclass
        end
        
        function ground = compute_diagnostic(ground, forcing) %function from base class fully overwritten
            ground = compute_diagnostic@GROUND_base_class(ground);
        end
        
        function ground = getThermalPropertiesInfiltration(ground, forcing)
            %[c_temp, k_temp, k_eff, lwc_temp]
            %------- unused grid cells --------------------------------------------
            c_temp(GRID.air.cT_domain) = PARA.constants.c_a;
            k_temp(GRID.air.cT_domain) = PARA.constants.k_a;
            lwc_temp(GRID.air.cT_domain) = 0;
            
%             %------- soil domain --------------------------------------------------
%             [c_temp(GRID.soil.cT_domain),...
%                 k_temp(GRID.soil.cT_domain),...
%                 lwc_temp(GRID.soil.cT_domain)] = readThermalParameters(T(GRID.soil.cT_domain), GRID, PARA);
%             
%             %adjust for the unfrozen part of the domain
%             c_temp(GRID.soil.cT_domain) = double(T(GRID.soil.cT_domain)<=0).*c_temp(GRID.soil.cT_domain) + double(T(GRID.soil.cT_domain)>0).* capacityUnfrozen(wc,GRID,PARA);
%             k_temp(GRID.soil.cT_domain) = double(T(GRID.soil.cT_domain)<=0).*k_temp(GRID.soil.cT_domain) + double(T(GRID.soil.cT_domain)>0).* conductivityUnfrozen(wc,GRID,PARA);
%             lwc_temp(GRID.soil.cT_domain) = double(T(GRID.soil.cT_domain)<=0).*lwc_temp(GRID.soil.cT_domain) + double(T(GRID.soil.cT_domain)>0).* wc;
%             
%             %-------- set higher conductivity for free water ----------------------
%             % now done in conductivityUnfrozen
%             
%             
%             %------- snow domain --------------------------------------------------
%             c_temp(GRID.snow.cT_domain) = cap_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
%                 GRID.snow.Snow_w(GRID.snow.cT_domain),...
%                 GRID.snow.Snow_a(GRID.snow.cT_domain),...
%                 PARA);
%             
%             k_temp(GRID.snow.cT_domain) = cond_snow(GRID.snow.Snow_i(GRID.snow.cT_domain),...
%                 GRID.snow.Snow_w(GRID.snow.cT_domain),...
%                 GRID.snow.Snow_a(GRID.snow.cT_domain));
%             
%             lwc_temp(GRID.snow.cT_domain) = GRID.snow.Snow_w(GRID.snow.cT_domain)./GRID.general.K_delta(GRID.snow.cT_domain);
%             
%             %------- interpolate conductivity to K-grid ---------------------------
%             k_eff(2:end-1) = GRID.general.K_delta(1:end-1)./(2.*GRID.general.cT_delta) .* (1./k_temp(1:end-1)).^2 ...
%                 + GRID.general.K_delta(2:end)  ./(2.*GRID.general.cT_delta) .* (1./k_temp(2:end)).^2;
%             
%             k_eff(2:end-1) = k_eff(2:end-1).^(-0.5);
%             
%             k_eff(1)     = k_temp(1);
%             k_eff(end)   = k_temp(end);
%             
%             %------ correct upper most value below air-domain ---------------------
%             k_eff(GRID.air.K_domain_lb+1) = k_temp(GRID.air.K_domain_lb+1);
        end
        
    end
end

