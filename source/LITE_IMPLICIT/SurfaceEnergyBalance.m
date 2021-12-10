% clear all;
% close all;
%% 0. set arbitray state and forcing values
% surface temperature TODO: replace with actual surface temperature
% STATE.T = 16 .* ones(10,1);

% forcing variables TODO: replace with actual forcing variables
% FORCING.z = 2.0;
% FORCING.wind = 5.04;
% FORCING.Tair = 15.5;
% FORCING.Sin = 186.1;
% FORCING.Lin = 351.0;
% FORCING.q = 0.0092;
% FORCING.p = 100208.2;

%

%% 0. set parameters 
% % surface characteristics TODO: update depending on surface state (soil, snow, water, ice,...)
% PARA.alpha = 0.2;                                                               % surface albedo [-]
% PARA.epsilon = 0.97;                                                            % surface emissivity [-]
% PARA.z0 = 1e-3;                                                                 % surface roughness length [m]
% PARA.rs = 50.0;                                                                 % surface resistance against evapotranspiration and sublimation [s/m]
% 
% % natural "constants"
% PARA.sigma = 5.6704e-8;                                                         % Stefan-Boltzmann constant
% PARA.kappa = 0.4;                                                               % von Kármán constant [-]
% PARA.gamma = 0.622;                                                             % Psychrometric constant [-]
% PARA.R_a = 287.058;                                                             % specific gas constant of air [J/(kg*K)]
% PARA.g = 9.81;                                                                  % gravitational acceleration [m/s^2]
% 
% % material properties (assumed to be constant)
% PARA.rho_a = 1.293;                                                             % density of air at standard pressure and 0°C [kg/m^3]
% PARA.c_a = 1005.7;                                                              % volumetric heat capacity of dry air at standard pressure and 0°C [J/(m^3*K)]
% 
% % parameters for Monin-Obukhov theory and stability functions according to Businger (1971) and Byun (1990)
% PARA.Pr_0 = 0.74;                                                               % turbulent Prandtl number
% PARA.beta_m = 4.7;
% PARA.beta_h = PARA.beta_m/PARA.Pr_0;
% PARA.gamma_m = 15.0;
% PARA.gamma_h = 9.0;


%% call seb function to update STATE
% STATE.Tsurf = 20; 
% GRID.k = 1 * ones(10,1);
% GRID.Delta = 0.01 * ones(10,1);
% 
% STATE = surfaceEnergyBalance_Byun( STATE, PARA, FORCING, GRID );
% fprintf("Before iteration: Tsurf=%0.2f, Qg=%0.2f\n", STATE.Tsurf,STATE.Qg );
% 
% 
% delta = 1;
% i=0;
% while abs(delta)>1e-3 && i<1000
%     i=i+1;
%     Tsurf_old = STATE.Tsurf;
%     STATE = surfaceEnergyBalance_Byun( STATE, PARA, FORCING, GRID );
%     delta = STATE.Tsurf - Tsurf_old;
% end
% 
% if i<1000
%     fprintf("After %d iterations: Tsurf=%0.2f, Qg=%0.2f\n", i, STATE.Tsurf,STATE.Qg );
% else
%     fprintf("Did not converge.\n");
% end

%%
function SEB = SurfaceEnergyBalance( T, Water, Porosity, Zp, dxp, kp, SEB, PARA, FORCING, tile_idx, surf_idx, t)
    albedo = PARA.albedo(tile_idx);
    emissivity = PARA.emissivity(tile_idx);

    % 1. calculate radiation budget
    % outgoing shortwave radiation as reflected
    SEB.Sout(tile_idx) = -albedo * FORCING.Sin(t);                                                       % Eq. (2) in Westermann et al. (2016)
    % outgoing longwave radiation composed of emitted and reflected radiation
    SEB.Lout(tile_idx) = -emissivity * PARA.sigma * (SEB.Tsurf(tile_idx)+273.15).^4 - (1 - emissivity) * FORCING.Lin(t);  % Eq. (3) in Westermann et al. (2016)
    % net radiation budget
    SEB.Qnet(tile_idx) = FORCING.Sin(t) + SEB.Sout(tile_idx) + FORCING.Lin(t) + SEB.Lout(tile_idx);

    % 2. calcuate turbulent heat flux budget
    % determine atmospheric stability conditions
    SEB.Lstar(tile_idx) = Lstar_Byun(SEB, PARA, FORCING, tile_idx,t);
    SEB.ustar(tile_idx) = ustar(SEB, PARA, FORCING, tile_idx,t);
    % sensible heat flux
    SEB.Qh(tile_idx) = Qh(SEB, PARA, FORCING, tile_idx,t);
    % latent heat flux
    SEB.Qe(tile_idx) = Qe(T, Water, Porosity, Zp, dxp, SEB, PARA, FORCING, tile_idx,t);

    % 3. determine ground heat flux as the residual of the radiative and turbulent fluxes
    SEB.Qg(tile_idx) = SEB.Qnet(tile_idx) - SEB.Qh(tile_idx) - SEB.Qe(tile_idx);                                    % essentially Eq. (1) in Westermann et al. (2016)
    
    % 4. tranlate ground heat flux into a new surface temperature
    SEB.Tsurf(tile_idx) = T(surf_idx,tile_idx) + SEB.Qg(tile_idx) * dxp(surf_idx,tile_idx) / 2. / kp(surf_idx,tile_idx); 
    
    if SEB.Tsurf(tile_idx)<-100 || SEB.Tsurf(tile_idx)>100 
        warning("unrealistic Tsurf");
    end
end

%% auxiliary functions    
% Lstar according to analytical/approximate solution for zeta by Byun 1990
% using bulk Richardson number
function Lstar = Lstar_Byun(SEB,PARA,FORCING,tile_idx,t)
    g = PARA.g;
    R_a = PARA.R_a;
    c_p = PARA.c_a / PARA.rho_a;                                            % specific heat capacity of air at constant pressure
    Tz = FORCING.Tair(t)+273.15;                                               % air temperature at height z over surface
    T0 = SEB.Tsurf(tile_idx)+273.15;                                                % surface temperature
    p = FORCING.p(t);                                                          % atmospheric pressure at surface (height z)
    p0 = FORCING.p(t);                                                         % normal pressure (for now assumed to be equal to p)
    uz = FORCING.wind(t);                                                      % wind speed at height z
    z = FORCING.z;                                                          % height z of wind forcing
    z0 = PARA.z0(tile_idx);                                                           % aerodynamic roughness length [m]
    Pr_0 = PARA.Pr_0;                                                       % turbulent Prandtl number
    gamma_h = PARA.gamma_h;
    gamma_m = PARA.gamma_m;
    beta_h = PARA.beta_h;
    beta_m = PARA.beta_m;

    Theta = Tz * (p0/p)^(R_a/c_p);                                          % potential temperature (for now identical to actual temperature)
    Theta0 = T0 * (p0/p)^(R_a/c_p);

    % calcuate bulk Richardson number
    Ri_b = g / Theta0 * (Theta - Theta0) * (z - z0) / uz^2;                 % bulk Richardson number, eq. (9) in Byun 1990

    % calulate ζ
    a = (z / (z-z0)) * log(z/z0);
    if Ri_b>0 % stable conditions
        b = 2 * beta_h * (beta_m * Ri_b - 1);
        c = -(2 * beta_h * Ri_b - 1) - (1 + (4 * (beta_h - beta_m) * Ri_b) / Pr_0 )^0.5;
        zeta = a / b * c;                                                   % eq. (19) in Byun 1990
    else % unstable conditions
        s_b = Ri_b / Pr_0;                                                  % eq. (30) in Byun 1990
        Q_b = 1/9 * ( 1/gamma_m^2 + 3 * gamma_h/gamma_m * s_b^2 );          % eq. (31) in Byun 1990
        P_b = 1/54 * ( -2/gamma_m^3 + 9/gamma_m * (3 - gamma_h/gamma_m) * s_b^2 );              % eq. (32) in Byun 1990
        Theta_b = acos( P_b / Q_b^(3/2) );                                      % eq. (33) in Byun 1990
        T_b = ( (P_b^2 - Q_b^3)^0.5 + abs(P_b) )^(1/3);                     % eq. (34) in Byun 1990
        if Q_b^3-P_b^2 >= 0
            zeta = a * (-2 * Q_b^0.5 * cos(Theta_b/3) + 1/(3*gamma_m) );                % eq. (28) in Byun 1990
        else
            zeta = a * (-(T_b + Q_b/T_b ) + 1/(3*gamma_m) );                        % eq. (29) in Byun 1990
        end
    end

    % calculate Lstar
    Lstar = z / zeta;

    % upper and lower limits for Lstar
    if abs(Lstar) < 1e-7
        Lstar = sign(Lstar) * 1e-6;
    end
    if abs(Lstar) > 1e+7
        Lstar = sign(Lstar) * 1e+6;
    end
end

function ustar = ustar(SEB,PARA,FORCING,tile_idx,t)
    kappa = PARA.kappa;
    uz = FORCING.wind(t);                                           
    z = FORCING.z;                                                     
    z0 = PARA.z0(tile_idx);                                               
    %Lstar = STATE.Lstar;
    % calculate ustar
    ustar = kappa * uz ./ (log(z / z0) - Psi_M(z / SEB.Lstar(tile_idx), z0 / SEB.Lstar(tile_idx), PARA));     % Eq. (7) in Westermann et al. (2016)
end

%Saturation pressure of water/ice according to the empirical August-Roche-Magnus formula
%Note: T is passed [K] and converted to [°C]
% Eq. (B3) in Westermann et al. (2016)
function estar = estar(T)
    if T > 0
        estar = 611.2 * exp(17.62 * (T - 273.15) / (243.12 + (T - 273.15))) ; 
    else
        estar = 611.2 * exp(22.46 * (T - 273.15) / (272.62 + (T - 273.15))) ;
    end
end

% sensible turbulent heat flux
function Qh = Qh(SEB, PARA, FORCING,tile_idx,t)
    kappa = PARA.kappa;
    R_a = PARA.R_a;
    Tz = FORCING.Tair(t)+273.15;                                               % air temperature
    T0 = SEB.Tsurf(tile_idx)+273.15;                                                % surface temperature
    c_p = PARA.c_a / PARA.rho_a;                                            % specific heat capacity of air at constant pressure
    z = FORCING.z;                                                          % height at which forcing data are provided
    Lstar = SEB.Lstar(tile_idx);
    ustar = SEB.ustar(tile_idx);
    p = FORCING.p(t);
    z0 = PARA.z0(tile_idx);
    rho_a = p / (Tz * R_a);                                                 % density of air at surface air temperature and surface pressure [kg/m^3]
    % calculate aerodynamic resistance
    raH = (kappa * ustar).^(-1) * (log(z / z0) - Psi_HW(z / Lstar, z0 / Lstar, PARA));    % Eq. (6) in Westermann et al. (2016)
    % calcuate sensible heat flux
    Qh = -rho_a * c_p * (Tz - T0) / raH;                                    % Eq. (4) in Westermann et al. (2016)
end

% latent turbulent heat flux
function Qe = Qe(T, Water, Porosity, Zp, dxp, SEB, PARA, FORCING,tile_idx,t)
    kappa = PARA.kappa;
    gamma = PARA.gamma;
    R_a = PARA.R_a;
    Tz = FORCING.Tair(t)+273.15;                                                            % air temperature at height z over surface
    T0 = SEB.Tsurf(tile_idx)+273.15;                                                             % surface temperature
    p = FORCING.p(t);                                                                % atmospheric pressure at surface
    qz = FORCING.q(t);                                                               % specific humidity at height h over surface
    z = FORCING.z;                                                                % height at which forcing data are provided
    rs = PARA.rs(tile_idx);                                                                 % surface resistance against evapotranspiration / sublimation [1/m]
    Lstar = SEB.Lstar(tile_idx);
    ustar = SEB.ustar(tile_idx);
    z0 = PARA.z0(tile_idx);                                                                 % aerodynamic roughness length [m]
    rho_a = p / (Tz * R_a);                                                       % density of air at surface air temperature and surface pressure [kg/m^3]
    q0 = gamma * estar(T0) / p;                                                   % saturation pressure of water/ice at the surface; Eq. (B1) in Westermann et al (2016)
    % calculate aerodynamic resistance
    raW = (kappa * ustar).^(-1) * (log(z / z0) - Psi_HW(z / Lstar, z0 / Lstar, PARA));  % aerodynamic resistance Eq. (6) in Westermann et al. (2016)
    % determine latent heat of phase change
    if T0<=273.15  
        L = 0;%1000 * (2834.1 - 0.29 * (T0-273.15) - 0.004 * (T0-273.15)^2);    % for now: set subliamation and deposition to zero to be consistent
    else
        L = 1000 * (2500.8 - 2.36 * (T0-273.15) + 0.0016 * (T0-273.15)^2 - 0.00006 * (T0-273.15)^3);
    end
    %calculate Q_E
    Qe_pot = -rho_a * L * (qz - q0) / (raW) ; 
    if qz > q0 % condensation / deposition
        Qe = Qe_pot;                             % Eq. (5) in Westermann et al. (2016) % condensation (no aerodynamics resistance)
    else
        % here, we distinguish depending on whether there is snow, water,
        % ice, or soil at the surface
        %if snowcover or lake or ice (no soil) 
        %if SnowDepth > 0
        %    Qe = Qe_pot * ( raW / (raW + rs) ) ;       % evaporation (account for surface resistance against evapotranspiration/sublimation)
        %else
            % water availability factor
        b = beta( T, Water, Porosity, Zp, dxp, PARA);
        if b > 0
            Qe = Qe_pot * b;
        else
            Qe = Qe_pot * ( raW / (raW + rs) ) ;     
        end
        %end
    end
end

function beta = beta(T, Water, Porosity, Zp, dxp, PARA)
    evaporationDepth = PARA.evaporationDepth;
    fieldCapacity = Porosity .* PARA.fieldCapacityFactor;
    
    % depth weighting
    depthWeighting = exp( - Zp / evaporationDepth );
    depthWeighting( depthWeighting > 1 ) = 0;
    depthWeighting( depthWeighting < 0.05 ) = 0;
    depthWeighting = depthWeighting .* dxp / sum(depthWeighting .* dxp);    % normaliz, taking into account grid cell thickness
    
    % temperature weighting
    tempWeighting = double(T>0);    % only consider unfrozen cells
    
    % hydraulic weighting
    hydraulicWeighting = double(Water>=fieldCapacity) + double(Water<fieldCapacity).*0.25.*(1-cos(pi().*Water./fieldCapacity)).^2;

    beta = sum( depthWeighting .* tempWeighting .* hydraulicWeighting );
end
                     
%Integrated stability function for heat/water transport; Busigner 1971              
function res = Psi_HW(zeta1,zeta2,PARA)
    gamma_h = PARA.gamma_h;
    beta_h = PARA.beta_h;
    if zeta1 <= 0 % neutral and unstable conditions (according to Businger, 1971)
        res = 2*log( ( (1-gamma_h*zeta1)^0.5 + 1 ) / ( (1-gamma_h*zeta2)^0.5 + 1 ) );      % eq. (15) in Byun 1990
    else     % stable stratification (according to Busigner, 1971)
        res = -beta_h*(zeta1-zeta2);                                                 % eq. (13) in Byun 1990
    end
end

%Integrated stability function for momentum transport; Busigner 1971
function res = Psi_M(zeta1,zeta2,PARA)
    gamma_m = PARA.gamma_m;
    beta_m = PARA.beta_m;
    if zeta1 <= 0 % neutral and unstable conditions (according to Businger, 1971)
        x=(1-gamma_m*zeta1)^0.25;
        x0=(1-gamma_m*zeta2)^0.25;
        res = 2*log( (1+x)/(1+x0) ) + log( (1+x^2) / (1+x0^2) ) - 2*atan(x) + 2*atan(x) ;    % eq. (14) in Byun 1990
    else     % stable stratification (according to Busigner, 1971)
        res = -beta_m*(zeta1-zeta2);                                        % eq. (12) in Byun 1990
    end
end      
