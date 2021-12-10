%========================================================================
% CryoGrid GROUND class GROUND_all_constant
% initialize as for GROUND_freeW_seb, and then do nothing
% S. Westermann, October 2020
%========================================================================

classdef GROUND_all_constant_BGC < GROUND_all_constant

    properties
        BGC
        IA_BGC
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

        
        function ground = provide_PARA(ground)
            
%             ground.PARA.albedo = [];  %surface albedo [-]
%             ground.PARA.epsilon = []; % surface emissivity [-]
%             ground.PARA.z0 = []; %roughness length [m]
%             
%             ground.PARA.rootDepth = []; %e-folding constant of transpiration reduction with depth [m]
%             ground.PARA.evaporationDepth = []; %e-folding constant of evaporation reduction reduction with depth [m]
%             ground.PARA.ratioET = []; %fraction of transpiration of total evapotranspiration [-]
%             ground.PARA.hydraulicConductivity = [];  %saturated hydraulic conductivity [m/sec]
%             
%             ground.PARA.dt_max = [];  %maximum possible timestep [sec]
%             ground.PARA.dE_max = [];  %maximum possible energy change per timestep [J/m3]
            ground = provide_PARA@GROUND_all_constant(ground);
            
            ground.PARA.BGC_CLASS = [];
            ground.PARA.target_grid_cell_size = 0.05;
        end
        
        function ground = provide_STATVAR(ground)
            
%             ground.STATVAR.upperPos = []; % upper surface elevation [m]
%             ground.STATVAR.lowerPos = []; % lower surface elevation [m]
%             ground.STATVAR.layerThick = []; % thickness of grid cells [m]
%             
%             ground.STATVAR.waterIce = []; % total volume of water plus ice in a grid cell [m3]
%             ground.STATVAR.mineral = []; % total volume of minerals [m3]
%             ground.STATVAR.organic = []; % total volume of organics [m3]
%             ground.STATVAR.energy = [];  % total internal energy [J]
%             ground.STATVAR.XwaterIce = [];
%             
%             ground.STATVAR.T = [];  % temperature [degree C]
%             ground.STATVAR.water = [];  % total volume of water [m3]
%             ground.STATVAR.ice = []; %total volume of ice [m3]
%             ground.STATVAR.air = [];  % total volume of air [m3] - NOT USED
%             ground.STATVAR.thermCond = []; %thermal conductivity [W/mK]
%             ground.STATVAR.hydraulicConductivity = []; % hydraulic conductivity [m/sec]
%             
%             ground.STATVAR.Lstar = []; %Obukhov length [m]
%             ground.STATVAR.Qh = []; %sensible heat flux [W/m2]
%             ground.STATVAR.Qe = []; % latent heat flux [W/m2]
%             
%             ground.STATVAR.field_capacity = []; %field capacity in fraction of the total volume [-]
%             ground.STATVAR.excessWater = 0;  %water volume overtopping first grid cell (i.e. surface water [m3]
            ground = provide_STATVAR@GROUND_all_constant(ground);
        end
        
        function ground = provide_CONST(ground)
            
%             ground.CONST.L_f = []; % volumetric latent heat of fusion, freezing
%             ground.CONST.c_w = []; % volumetric heat capacity water
%             ground.CONST.c_i = []; % volumetric heat capacity ice
%             ground.CONST.c_o = []; % volumetric heat capacity organic
%             ground.CONST.c_m = []; % volumetric heat capacity mineral
%             
%             ground.CONST.k_a = [];   % thermal conductivity air
%             ground.CONST.k_w = [];   % thermal conductivity water
%             ground.CONST.k_i = [];   % thermal conductivity ice
%             ground.CONST.k_o = [];   % thermal conductivity organic
%             ground.CONST.k_m = [];   % thermal conductivity mineral
%             
%             ground.CONST.sigma = []; %Stefan-Boltzmann constant
%             ground.CONST.kappa = []; % von Karman constant
%             ground.CONST.L_s = [];  %latent heat of sublimation, latent heat of evaporation handled in a dedicated function
%             
%             ground.CONST.cp = []; %specific heat capacity at constant pressure of air
%             ground.CONST.g = []; % gravitational acceleration Earth surface
%             
%             ground.CONST.rho_w = []; % water density
%             ground.CONST.rho_i = []; %ice density

              ground = provide_CONST@GROUND_all_constant(ground);
              
        end
        
        function ground = finalize_init(ground, tile) 

            ground = finalize_init@GROUND_all_constant(ground, tile);
            %ground.STATVAR.XwaterIce = ground.STATVAR.waterIce .* 0; 
            ground.STATVAR.T = ground.STATVAR.waterIce .* 0+5;
            
            class_handle = str2func(ground.PARA.BGC_CLASS);
            %ground.BGC = class_handle(-1,0,0,0); 
            ground.BGC = class_handle(); 
            ground.BGC.PARENT = ground;
            %remove this in the end
            ground.BGC = provide_PARA(ground.BGC);
            ground.BGC = provide_STATVAR(ground.BGC);
            ground.BGC = provide_CONST(ground.BGC);
            ground.BGC = finalize_init(ground.BGC, tile);
            
            ground.IA_BGC = IA_BGC_simple();
            ground.IA_BGC.BGC = ground.BGC;
            ground.IA_BGC.GROUND = ground;
            ground.BGC.IA_BGC = ground.IA_BGC;
            finalize_init(ground.IA_BGC, tile);

        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            
            ground.BGC = get_boundary_condition_u(ground.BGC, tile);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            
           ground.BGC = get_boundary_condition_l(ground.BGC, tile); 
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
           ground.BGC = get_derivatives_prognostic(ground.BGC, tile); 
        end
        
        function timestep = get_timestep(ground, tile) 
            timestep_ground = get_timestep@GROUND_all_constant(ground, tile);
            timestep_BGC = get_timestep(ground.BGC, tile);
            timestep = min(timestep_ground, timestep_BGC);
        end
        
        function ground = advance_prognostic(ground, tile)    
            
            ground.BGC = advance_prognostic(ground.BGC, tile);
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            
            ground.BGC = compute_diagnostic_first_cell(ground.BGC, tile);
        end
       
        function ground = compute_diagnostic(ground, tile)
            
            ground.BGC = compute_diagnostic(ground.BGC, tile);
                        
            ground = compute_diagnostic@GROUND_all_constant(ground, tile);
            
%             ground.STATVAR.water(1) = ground.STATVAR.layerThick(1).* ground.STATVAR.area(1)-ground.STATVAR.mineral(1)-ground.STATVAR.organic(1);
%             ground.STATVAR.waterIce(1) = ground.STATVAR.layerThick(1).* ground.STATVAR.area(1)-ground.STATVAR.mineral(1)-ground.STATVAR.organic(1);
        end
        
        function ground = check_trigger(ground, tile)
            
            ground.BGC = check_trigger(ground.BGC, tile);
            if size(ground.STATVAR.water,1) <size(ground.STATVAR.waterIce,1)
                ground.STATVAR.water = [ground.STATVAR.waterIce(1); ground.STATVAR.water];
                ground.STATVAR.T = [ground.STATVAR.T(1); ground.STATVAR.T];
                ground.STATVAR.mineral = [0; ground.STATVAR.mineral];
                ground.STATVAR.non_BGC_layerThick = [0; ground.STATVAR.non_BGC_layerThick]; %layerThick of each grid celll before BGC accumulation
                ground.STATVAR.non_BGC_mineral = [0; ground.STATVAR.non_BGC_mineral]; %layerThick of each grid celll before BGC accumulation
                ground.STATVAR.non_BGC_organic = [0; ground.STATVAR.non_BGC_organic];
            end
            if size(ground.STATVAR.water,1) >size(ground.STATVAR.waterIce,1)
                ground.STATVAR.water(1,:) = [];
                ground.STATVAR.T(1,:) = [];
                ground.STATVAR.mineral(1,:) = [];
                ground.STATVAR.non_BGC_layerThick(1,:) = []; %layerThick of each grid celll before BGC accumulation
                ground.STATVAR.non_BGC_mineral(1,:) = []; %layerThick of each grid celll before BGC accumulation
                ground.STATVAR.non_BGC_organic(1,:) = [];
            end
        end
        


    end
    
end
