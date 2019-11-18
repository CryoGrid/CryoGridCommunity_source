% inherits from GROUND_base_class and adds SEB as upper boundary
% no interaction with snow is possible here

classdef GROUND_vegetation_snow < GROUND_vegetation
    properties
        IA_CHILD
    end
    
    % Here comes the documentation of your new module
    % classdef GROUND_vegetation (Simone Stünzi, 20191010)
    % properties
    % CONST %constants
    % PARA %external service parameters, all other
    % STATVAR %state variables - choose wisely
    % TEMP %derivatives in prognostic
    % PREVIOUS %pointer to previous module
    % NEXT %pointer to next module
    % IA_PREVIOUS %pointer to interaction with previous module
    % IA_NEXT %pointer to interaction with next module
    %%%%%
    
    methods
        
        %mandatory functions for each class
        
        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            ground = provide_variables@GROUND_vegetation(ground); %call function of the base class
            ground = provide_PARA(ground); %add additional variables
            ground = provide_CONST(ground);
            ground = provide_STATVAR(ground);
        end
        
        function ground = assign_global_variables(ground, forcing)
            ground = assign_global_variables@GROUND_vegetation(ground, forcing); %call function of the base class
            ground.PARA.airT_height = forcing.PARA.airT_height;
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            ground = initialize_STATVAR_from_file@GROUND_vegetation(ground, grid, forcing, depths);
            ground = set_up_canopy(ground);
            ground = finalize_STATVAR(ground, forcing); %assign all variables, that must be calculated or assigned otherwise
        end
        
        function [ground] = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = get_boundary_condition_u(ground.IA_CHILD, forcing); %call boundary condition for child
            else
                ground = surface_energy_forest(ground, forcing);
            end
            %
            
            
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
            %ground = get_boundary_condition_l@GROUND_base_class(ground, forcing);
        end
        
        
        function ground = get_derivatives_prognostic(ground)
            %ground = get_derivatives_prognostic@GROUND_base_class(ground);
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = get_derivative_energy(ground.IA_CHILD); % boundary condition for child, omacts ground.TEMP.F_ub
            end
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
            %timestep = get_timestep@GROUND_base_class(ground);
            timestep = 24*60*60;
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = advance_prognostic(ground.IA_CHILD, timestep); %call function for child
            end
            
            %ground = advance_prognostic@GROUND_base_class(ground, timestep);
            %             ground.STATVAR.current_t = forcing.TEMP.t + timestep;
        end
        
                function ground = compute_diagnostic_first_cell(ground, forcing)
                    %ground = L_star(ground, forcing);
                end
        
        function ground = compute_diagnostic(ground, forcing)
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = compute_diagnostic(ground.IA_CHILD); %call function for child
            end
            
            %             ground = compute_diagnostic@GROUND_base_class(ground, forcing);
        end
        

        
        function ground = surface_energy_forest(ground, forcing)
            if forcing.TEMP.t >= ground.STATVAR.execution_t
                datestr(forcing.TEMP.t)
                
                vegetation = ground.STATVAR.vegetation;
                
                [vegetation] = set_up_forcing(vegetation, forcing);
%                 
%                 disp('Albedo used for canopy fluxes multilayer:')
%                 vegetation.flux.albsoib
                
                [vegetation] = canopy_fluxes_multilayer(vegetation);
                                
%                 [vegetation] = figures(vegetation);
                
                vegetation.flux.albsoib
                
                ground.STATVAR.vegetation = vegetation;
                
                ground.STATVAR.execution_t = ground.STATVAR.execution_t + 1/24;
                
                
                ground.TEMP.F_ub = ground.STATVAR.vegetation.mlcanopyinst.gsoi;
%                 ground.STATVAR.T(end) = ground.STATVAR.vegetation.mlcanopyinst.tveg(:,end,1)-273.15;
%                 ground.STATVAR.thermCond(end) = 0.025; % Thermal conductivity of air
%                 ground.STATVAR.layerThick(end) = ground.STATVAR.vegetation.mlcanopyinst.tveg(:,2);
%                 ground.STATVAR.Lout = ground.STATVAR.vegetation.mlcanopyinst.ircan;
%                 ground.STATVAR.Qh = ground.STATVAR.vegetation.mlcanopyinst.shveg;
%                 ground.STATVAR.Qe = ground.STATVAR.vegetation.mlcanopyinst.lhveg;
                
            end
            
        end
        
    end
    
end
