% inhertits from SNOW_base_class and adds SEB as upper boundary
% designed to function as CHILD of a GROUND class that is compatible with
% SNOW classes; compatible with interaction classes IA_SNOW_GROUND and IA_SNOW_GROUND_fcSimple_salt

classdef SNOW_seb_simple_vegetation_below < GROUND_base_class
    
    methods
        
        %mandatory functions for each class
        
        function snow = provide_variables(ground)  %initializes the subvariables as empty arrays
            snow = provide_variables@GROUND_base_class(ground); %call function of the base class
            snow = provide_PARA(ground); %add additional variables
            snow = provide_CONST(ground);
            snow = provide_STATVAR(ground);
        end
        
        function snow = assign_global_variables(ground, forcing)
            snow = assign_global_variables@GROUND_base_class(ground, forcing); %call function of the base class
            snow.PARA.airT_height = forcing.PARA.airT_height;
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            ground = initialize_STATVAR_from_file@GROUND_base_class(ground, grid, forcing, depths);
            ground = set_up_canopy(ground);
            ground = finalize_STATVAR(ground, forcing); %assign all variables, that must be calculated or assigned otherwise
        end
        
        function [ground] = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            if forcing.TEMP.t >= ground.STATVAR.execution_t
                datestr(forcing.TEMP.t)
                
                vegetation = ground.STATVAR.vegetation;
                
                vegetation.flux.albsoib = [0.8,0.8]; % Direct beam albedo of ground (soil)
                vegetation.flux.albsoid = [0.8,0.8]; % Diffuse albedo of ground (soil)
                
                [vegetation] = set_up_forcing(vegetation, forcing);
                
                [vegetation] = canopy_fluxes_multilayer(vegetation);
                
                %                                 [vegetation] = figures(vegetation);
                
                vegetation.flux.albsoib
                
                ground.STATVAR.vegetation = vegetation;
                
                ground.STATVAR.execution_t = ground.STATVAR.execution_t + 1/24;
                
                
                ground.TEMP.F_ub = ground.STATVAR.vegetation.mlcanopyinst.gsoi;
                ground.STATVAR.T(end) = ground.STATVAR.vegetation.mlcanopyinst.tveg(:,2,1)-273.15;
                ground.STATVAR.thermCond(end) = 0.025; % Thermal conductivity of air
                ground.STATVAR.layerThick(end) = ground.STATVAR.vegetation.mlcanopyinst.tveg(:,2);
                ground.STATVAR.Lout = ground.STATVAR.vegetation.mlcanopyinst.ircan;
                ground.STATVAR.Qh = ground.STATVAR.vegetation.mlcanopyinst.shveg;
                ground.STATVAR.Qe = ground.STATVAR.vegetation.mlcanopyinst.lhveg;
                
            end
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = get_boundary_condition_u(ground.IA_CHILD, forcing); %call boundary condition for child
            end
            
            
            
        end
        
        %         function ground = get_boundary_condition_l(ground, forcing)
        %             %ground = get_boundary_condition_l@GROUND_base_class(ground, forcing);
        %           end
        
        
        function ground = get_derivatives_prognostic(ground)
            %ground = get_derivatives_prognostic@GROUND_base_class(ground);
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = get_derivative_energy(ground.IA_CHILD); % boundary condition for child, omacts ground.TEMP.F_ub
            end
        end
        
        
        
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
            %timestep = get_timestep@GROUND_base_class(ground);
            timestep = 24*60*60;
            ground.STATVAR.vegetation.flux.albsoib
        end
        
        function ground = advance_prognostic(ground, timestep) %real timestep derived as minimum of several classes in [sec] here!
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = advance_prognostic(ground.IA_CHILD, timestep); %call function for child
            end
            
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
        % % %         function ground = surface_energy_balance(ground, forcing)
        % % %             % Set up required forcing structure
        % % %             [ground] = SetUpForcing(ground, forcing);
        % % %             % Run the Bonan Model
        % % %             [ground] = CanopyFluxesMultilayer(ground); % , counter
        % % %             % Si:  ground.TEMP.F_ub = ground.STATVAR.vegetation.[GroundHeatFlux (da wo er jetzt halt grad drinn ist in meinem Bonan Vegi Modell)]
        % % %             ground.TEMP.F_ub = ground.STATVAR.vegetation.mlcanopyinst.gsoi; %ground.STATVAR.vegetation.mlcanopyinst.gsoi;
        % % % %             % Si:  ground.TEMP.F_ub = ground.STATVAR.vegetation.[GroundHeatFlux (da wo er jetzt halt grad drinn ist in meinem Bonan Vegi Modell)]
        % % % %             % ground.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - ground.STATVAR.Lout - ground.STATVAR.Sout - ground.STATVAR.Qh - ground.STATVAR.Qe;
        % % %         end
        
        
        
        function snow = check_trigger(snow)
            if size(snow.STATVAR.energy,1) ==1 && snow.STATVAR.ice < snow.PARA.swe_per_cell./2
                ground = snow.NEXT;
                
                ground.PREVIOUS = snow.PREVIOUS; %reassign ground
                ground.PREVIOUS.NEXT = ground;
                ground.IA_PREVIOUS=[];
                
                snow.NEXT =[]; %reassign snow -> cut all dependencies
                snow.PREVIOUS =[];
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