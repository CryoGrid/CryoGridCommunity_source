%Sediment with heat and salt diffusion

classdef GROUND_Sediment_heat_fcurve_saltDiff_sedimentation < GROUND_Sediment_heat_fcurve_saltDiff
    properties
        IA_CHILD
        IA_PARENT
    end
    methods
        
        %mandatory functions for each class
        
        function xls_out = write_excel(ground)
            xls_out = {'CLASS','index',NaN,NaN,NaN;'GROUND_Sediment_heat_fcurve_saltDiff_sedimentation',1,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN;NaN,'value','default','unit',NaN;'albedo',0.200000000000000,0.150000000000000,'[-]','surface albedo';'epsilon',0.990000000000000,0.990000000000000,'[-]','surface emissivity';'z0',0.00100000000000000,0.00100000000000000,'[m]','roughness length';'rootDepth',0.200000000000000,NaN,'[m]',NaN;'evaporationDepth',0.0500000000000000,NaN,'[m]',NaN;'ratioET',0.500000000000000,NaN,'[-]',NaN;'hydraulicConductivity',1.00000000000000e-05,NaN,'[m/sec]',NaN;'    ','    ','    ',NaN,'    ';'dt_max',3600,3600,'[sec]','longest possible timestep';'dE_max',50000,50000,'[J/m3]','maximum change of energy per timestep';'CLASS_END',NaN,NaN,NaN,NaN};
        end

        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            ground = provide_variables@GROUND_Sediment_heat_fcurve_saltDiff(ground); %call function of the base class
        end
        
        function ground = assign_global_variables(ground, forcing)
            ground = assign_global_variables@GROUND_Sediment_heat_fcurve_saltDiff(ground, forcing); %call function of the base class
        end
        
        function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
            %move upperPos down to account for sedimentation
            totalSedimentation = 20e-5*abs(forcing.DATA.timeForcing(1));
            ground.STATVAR.upperPos = ground.STATVAR.upperPos - totalSedimentation;
            
            ground = initialize_STATVAR_from_file@GROUND_Sediment_heat_fcurve_saltDiff(ground, grid, forcing, depths);
        end
        
        function ground = initialize_zero_depth(ground, parentGround)
            %remove all pointer
            ground.IA_NEXT = [];
            ground.IA_CHILD = [];
            ground.NEXT = [];
            ground.PREVIOUS = [];
            
            %this function is here to overwrite everything with zeros that
            %is already initialized when copying a new class. 
            ground.STATVAR.T = parentGround.STATVAR.T(1);
            ground.STATVAR.c_eff = parentGround.STATVAR.c_eff(1);
            ground.STATVAR.saltConc = parentGround.STATVAR.saltConc(1);
            ground.STATVAR.porosity = parentGround.STATVAR.porosity(1);

            ground.STATVAR.water = parentGround.STATVAR.water(1); % [m]
            ground.STATVAR.mineral = parentGround.STATVAR.mineral(1); % [m]
            ground.STATVAR.organic = parentGround.STATVAR.organic(1); % [m]
            ground.STATVAR.ice = parentGround.STATVAR.ice(1);
            ground.STATVAR.air = parentGround.STATVAR.air(1);  % [m]
            ground.STATVAR.thermCond = parentGround.STATVAR.thermCond(1);
            
            ground.STATVAR.liqWater = parentGround.STATVAR.liqWater(1); % [m]
            ground.STATVAR.saltConc = ground.STATVAR.liqWater*ground.PARA.salinity;
            
            ground.STATVAR.waterIce = parentGround.STATVAR.waterIce_tot(1);
            ground.STATVAR.energy = parentGround.STATVAR.energy(1);
            ground.STATVAR.deltaT = parentGround.STATVAR.deltaT(1);
            ground.STATVAR.waterIce_tot = parentGround.STATVAR.waterIce_tot(1);
            ground.STATVAR.mineral_tot = parentGround.STATVAR.mineral_tot(1);
            ground.STATVAR.organic_tot = parentGround.STATVAR.organic_tot(1);
            ground.STATVAR.energy_tot = parentGround.STATVAR.energy_tot(1);
            ground.STATVAR.water_tot = parentGround.STATVAR.water_tot(1);
            ground.STATVAR.ice_tot = parentGround.STATVAR.ice_tot(1);
            ground.STATVAR.air_tot = parentGround.STATVAR.air_tot(1);
                        
            ground.STATVAR.Tmelt = parentGround.STATVAR.Tmelt(1);
            ground.STATVAR.saltDiff = parentGround.STATVAR.saltDiff(1);
            
            ground.STATVAR.layerThick = 0; % [m]
            ground.STATVAR.upperPos = parentGround.STATVAR.upperPos;
            ground.STATVAR.lowerPos = parentGround.STATVAR.upperPos;
            
            ground.TEMP.sedimentation = 0;
        end
        
        
        function ground = get_boundary_condition_u(ground, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
        
            ground = get_boundary_condition_u@GROUND_Sediment_heat_fcurve_saltDiff(ground, forcing);       
            
            %child boundary condition gives a new boundary condition
            %(similar to IA_HEAT/_SALT) to the parent if STATUS == 1
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = get_boundary_condition_u(ground.IA_CHILD, forcing); %call boundary condition for child
            end
                  
        end
        
        function ground = get_boundary_condition_l(ground, forcing)
            ground = get_boundary_condition_l@GROUND_Sediment_heat_fcurve_saltDiff(ground, forcing);
        end
        
        
        function ground = get_derivatives_prognostic(ground)
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = get_derivative_temperature_salt(ground.IA_CHILD); % boundary condition for child, omacts ground.TEMP.F_ub
            end
            
            ground = get_derivatives_prognostic@GROUND_Sediment_heat_fcurve_saltDiff(ground);
                       
        end
        
        function timestep = get_timestep(ground)  %could involve check for several state variables
           timestep = get_timestep@GROUND_Sediment_heat_fcurve_saltDiff(ground);
           if ~isempty(ground.IA_CHILD)
                timestep_child = get_timestep@GROUND_Sediment_heat_fcurve_saltDiff(ground.IA_CHILD.IA_CHILD);
                timestep = min(timestep, timestep_child);
                timestep = max(timestep, 0.0001);
           end
            
        end
        
        function ground = advance_prognostic(ground, timestep) 
            ground = advance_prognostic@GROUND_Sediment_heat_fcurve_saltDiff(ground, timestep); %advance energy and route down water
            
            if ~isempty(ground.IA_CHILD)
                ground.IA_CHILD = advance_prognostic(ground.IA_CHILD, timestep); %call function for child CHECK!!
            end
        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)
            ground = compute_diagnostic_first_cell@GROUND_Sediment_heat_fcurve_saltDiff(ground, forcing);
            %add birthing of new child here 
            if isempty(ground.IA_CHILD) && forcing.TEMP.surfaceState ~= 2
               %make new empyt child with status=0
               ground.IA_CHILD = IA_SEDIMENTATION();
               
               CURRENT = ground.IA_CHILD; %change to interaction class
               CURRENT.STATUS = 0; %sedimentation initially inactive
               CURRENT.IA_PARENT = ground;
               CURRENT.IA_CHILD = copy(ground); %child is a copy of current sediment class
               CURRENT.IA_CHILD.IA_PARENT = CURRENT;
                           
               %set parameters according to surfaceState
               switch forcing.TEMP.surfaceState
                   case 0 %subsea
                       CURRENT.IA_CHILD.STATVAR.soilType = 0;
                       CURRENT.IA_CHILD.PARA.salinity = 895;
                       CURRENT.IA_CHILD.PARA.alpha = 6.50e-01;
                       CURRENT.IA_CHILD.PARA.n = 1.67;
                       CURRENT.IA_CHILD.PARA.delta = 1.00e-07;
                       CURRENT.IA_CHILD.PARA.saltDiff0 = 8.00e-10;
                       CURRENT.IA_CHILD.PARA.sedRate = 30e-5;
                                  
                       CURRENT.IA_CHILD.TEMP.surfaceState = forcing.TEMP.surfaceState;

                   case 1 %subaerial
                       CURRENT.IA_CHILD.STATVAR.soilType = 1;
                       CURRENT.IA_CHILD.PARA.salinity = 0;
                       CURRENT.IA_CHILD.PARA.alpha = 4.06;
                       CURRENT.IA_CHILD.PARA.n = 2.03;
                       CURRENT.IA_CHILD.PARA.delta = 1.00e-07;
                       CURRENT.IA_CHILD.PARA.saltDiff0 = 8.00e-10;
                       CURRENT.IA_CHILD.PARA.sedRate = 10e-5;
                       
                       CURRENT.IA_CHILD.TEMP.surfaceState = forcing.TEMP.surfaceState;
               end
               
               %overwrite everything in the copied class
               %most statevariables get inherited from the parent, salinity
               %comes from the paramters
               CURRENT.IA_CHILD = initialize_zero_depth(CURRENT.IA_CHILD, CURRENT.IA_PARENT);
               ground.IA_CHILD = CURRENT; %reassign to ground
               
            end
        end
        
        function ground = compute_diagnostic(ground, forcing)
            if ~isempty(ground.IA_CHILD)
            	%if child exceeds threshold, make child a new module,
                %move to top, and update forcing.TEMP.elev
                %add isostasy here, maybe work with flag to trigger isostatic
                %movement in the lower modules
                ground.IA_CHILD = compute_diagnostic(ground.IA_CHILD); %call function for interaction
                ground.IA_CHILD = check_trigger(ground.IA_CHILD);
                
                if ~isempty(ground.IA_PREVIOUS) %if the child is the new upper module,     
                                             %the previous pointer is no
                                             %longer empty
                    forcing.TEMP.elev = ground.PREVIOUS.STATVAR.upperPos;
                end
            end
            
            ground = compute_diagnostic@GROUND_Sediment_heat_fcurve_saltDiff(ground, forcing);
        end
        
        function ground = troubleshoot(ground)
            ground = checkNaN(ground);
        end
        
    
    end
    
end
