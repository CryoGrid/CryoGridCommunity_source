function ground = provide_CONST(ground)
    
    ground.CONST.L_f = [];
    
    ground.CONST.c_w = [];
    ground.CONST.c_i = []; 
    ground.CONST.c_o = []; 
    ground.CONST.c_m = []; 
    
    ground.CONST.k_a = [];       %air [Hillel(1982)]
    ground.CONST.k_w = [];        %water [Hillel(1982)]
    ground.CONST.k_i = [];         %ice [Hillel(1982)]
    ground.CONST.k_o = [];        %organic [Hillel(1982)]
    ground.CONST.k_m = [];
   
    
    ground.CONST.kappa = 0.4;                                     % von Kármán constant [-]
    ground.CONST.sigma = 5.6704e-8;                               % Stefan-Boltzmann constant [ W / (m^2 K^4) ]
    ground.CONST.g = 9.81;                                        % gravitational acceleration [m/s^2]
    ground.CONST.p_0 = 100500;                                    % normal pressure (sea level) [Pa=kg/(m s^2)]
    %water
    ground.CONST.rho_w = 1000;                                    % density of liquid water (and ice) [kg/m^3]
    ground.CONST.c_w = 4200 * ground.CONST.rho_w;               % volumetric heat capacity of water [ J / (m^3 K) ]
    ground.CONST.k_w = 0.57;                                      % heat conductivity of water [ W/(mK) ] [Hillel(1982)]
    %ice
    ground.CONST.rho_i = 1000;                                    % density of ice, assumed to be equal to that of water [kg/m^3]
    ground.CONST.c_i = 1900 * ground.CONST.rho_i;               % volumetric heat capacity of ice [ J / (m^3 K) ]
    ground.CONST.k_i = 2.2;                                       % heat conductivity of ice [ W/(mK) ] [Hillel(1982)]
    %latent heat of water
    ground.CONST.T_f = 273.15;                                    % freezing point of water / zero degree Celsius [K]
    ground.CONST.L_sl = 334e3;                                    % specific latent heat of fusion of water [J/kg]            [AMS]
    ground.CONST.L_lg = 2501e3;                                   % specific latent heat of vaporization of water [J/kg]      [AMS]
    ground.CONST.L_sg = ground.CONST.L_sl + ground.CONST.L_lg;% specific latent heat of sublimation of water [J/kg]       
    %air
    ground.CONST.rho_a = 1.293;                                   % density of air [kg/m^3] @ 0°C, sea level
    ground.CONST.c_a = 1005.7 * ground.CONST.rho_a;             %c_a= 0.00125*10^6;%[J/m^3 K]   % volumetric heat capacity of dry air [J/(m^3 K)] @ 0°C, sea level, isobar
    ground.CONST.k_a = 0.0243;                                    %ka=0.025; [Hillel(1982)]       % heat conductivity of air [ W/(mK)] @ 0 °C, sea level pressure     
    ground.CONST.R_a = 287.058;                                   % specific gas constant of air [ J/(kg K) ]
    % organic
    %ground.CONST.rho_o = 1; % n.a.
    ground.CONST.c_o = 2.5e6; %[J/(K m^3)]                        % volumetric heat capacity of organic material [J/(K m^3)]
    ground.CONST.k_o = 0.25;                                      % heat conductivity of organic material [ W/(mK) ] [Hillel(1982)]
    % mineral
    %ground.CONST.rho_m = 1; % n.a.
    ground.CONST.c_m = 2.0e6; %[J/(K m^3)]                        % volumetric heat capacity of minearal material [J/(K m^3)]  
   ground.CONST.kh_bedrock =3;
   ground.CONST.k_m = ground.CONST.kh_bedrock;       

   
    