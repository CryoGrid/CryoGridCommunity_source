function [vegetation] = NormanRadiation (rho, tau, omega, td, tb, tbj, vegetation) %, p, ic
%   rho, tau, omega, tb, td, tbj, surfalb_inst, mlcanopy_inst)
%     % %DESCRIPTION:
%     % Compute solar radiation transfer through canopy using Norman (1979)

%     % %ARGUMENTS:
%     implicit none
%     type(bounds_type), intent(in) :: bounds
%     integer,  intent(in) :: num_exposedvegp                          % Number of non-snow-covered veg points in CLM patch filter
%     integer,  intent(in) :: filter_exposedvegp(:)                    % CLM patch filter for non-snow-covered vegetation
%     real(r8), intent(in) :: rho(bounds%begp:bounds%endp,1:numrad)    % Leaf/stem reflectance
%     real(r8), intent(in) :: tau(bounds%begp:bounds%endp,1:numrad)    % Leaf/stem transmittance
%     real(r8), intent(in) :: omega(bounds%begp:bounds%endp,1:numrad)  % Leaf/stem scattering coefficient
%     real(r8), intent(in) :: tb(bounds%begp:bounds%endp,1:nlevcan)    % Exponential transmittance of direct beam through a single leaf layer
%     real(r8), intent(in) :: td(bounds%begp:bounds%endp,1:nlevcan)    % Exponential transmittance of diffuse through a single leaf layer
%     real(r8), intent(in) :: tbj(bounds%begp:bounds%endp,0:nlevcan)   % Exponential transmittance of direct beam onto canopy layer
%     type(surfalb_type), intent(in) :: surfalb_inst
%     type(mlcanopy_type), intent(inout) :: mlcanopy_inst
%     %
%     % %LOCAL VARIABLES:
%     integer  :: f                                                % Filter index
%     integer  :: p                                                % Patch index for CLM g/l/c/p hierarchy
%     integer  :: c                                                % Column index for CLM g/l/c/p hierarchy
%     integer  :: ic                                               % Aboveground layer index
%     integer  :: icm1                                             % Layer below ic (ic-1)
%     integer  :: ib                                               % Waveband index
%     real(r8) :: suminc                                           % Incident radiation for energy conservation check
%     real(r8) :: sumref                                           % Reflected radiation for energy conservation check
%     real(r8) :: sumabs                                           % Absorbed radiation for energy conservation check
%     real(r8) :: error                                            % Error check
%     real(r8) :: trand                                            % Term for diffuse radiation transmitted by layer
%     real(r8) :: refld                                            % Term for diffuse radiation reflected by layer
%     integer  :: m                                                % Index to the tridiagonal matrix
%     real(r8) :: aic, bic                                         % Intermediate terms for tridiagonal matrix
%     real(r8) :: eic, fic                                         % Intermediate terms for tridiagonal matrix
neq = (vegetation.mlcanopyinst.ncan+1)*2;                  % Number of tridiagonal equations to solve
a = zeros(1,neq);
b = zeros(1,neq);                     % Entries in tridiagonal matrix
c = zeros(1,neq);
d = zeros(1,neq);                     % Entries in tridiagonal matrix
u = zeros(1,neq);                                % Tridiagonal solution

swup = zeros(1,vegetation.mlcanopyinst.ncan,vegetation.mlcanopyinst.nleaf); %(bounds%begp:bounds%endp,0:nlevcan,numrad)   % Upward diffuse solar flux above canopy layer (W/m2 ground)
swdn = zeros(1,vegetation.mlcanopyinst.ncan,vegetation.mlcanopyinst.nleaf); %(bounds%begp:bounds%endp,0:nlevcan,numrad)   % Downward diffuse solar flux onto canopy layer (W/m2 ground)

% Compute solar radiation transfer through canopy using Norman (1979)

% -----------------------------------------------------------------------
% Input
% rho            % Leaf reflectance
% tau            % Leaf transmittance
% omega          % Leaf scattering coefficient
% td             % Exponential transmittance of diffuse radiation through a single leaf layer
% tb             % Exponential transmittance of direct beam radiation through a single leaf layer
% tbcum          % Cumulatice exponential transmittance of direct beam onto a canopy layer

% vegetation.params.numrad = vegetation.params.numrad;  % Number of wavebands
% vegetation.params.npts   = vegetation.params.npts;    % Number of grid points to process
% vegetation.params.sun   = vegetation.params.sun  ;   % Index for sunlit leaf
% vegetation.params.sha   = vegetation.params.sha  ;   % Index for shaded leaf
% vegetation.canopy.ntop   = vegetation.canopy.ntop  ;  % Index for top leaf layer
% vegetation.canopy.nbot   = vegetation.canopy.nbot ;   % Index for bottom leaf layer
% vegetation.canopy.nsoi   = vegetation.canopy.nsoi  ;  % First canopy layer is soil
% vegetation.canopy.dlai   = vegetation.canopy.dlai  ;  % Layer leaf area index (m2/m2)
% vegetation.atmos.swskyb   = vegetation.atmos.swskyb ;  % Atmospheric direct beam solar radiation (W/m2)
% vegetation.atmos.swskyd   = vegetation.atmos.swskyd ;  % Atmospheric diffuse solar radiation (W/m2)
% vegetation.flux.fracsun   = vegetation.flux.fracsun ;  % Sunlit fraction of canopy layer
% vegetation.flux.fracsha   = vegetation.flux.fracsha ;  % Shaded fraction of canopy layer
% vegetation.flux.albsoib   = vegetation.flux.albsoib ;  % Direct beam albedo of ground (soil)
% vegetation.flux.albsoid   = vegetation.flux.albsoid ;  % Diffuse albedo of ground (soil)

% Output
% vegetation.flux.swleaf   = vegetation.flux.swleaf ;   % Leaf absorbed solar radiation (W/m2 leaf)
% vegetation.flux.swveg   = vegetation.flux.swveg   ;  % Absorbed solar radiation, vegetation (W/m2)
% vegetation.flux.swvegsun   = vegetation.flux.swvegsun ; % Absorbed solar radiation, sunlit canopy (W/m2)
% vegetation.flux.swvegsha   = vegetation.flux.swvegsha ; % Absorbed solar radiation, shaded canopy (W/m2)
% vegetation.flux.swsoi   = vegetation.flux.swsoi   ;  % Absorbed solar radiation, ground (W/m2)
% vegetation.flux.albcan   = vegetation.flux.albcan ;   % Albedo above canopy

% -----------------------------------------------------------------------

isun = vegetation.params.sun; % Array index for sunlit leaf
isha = vegetation.params.sha; % Array index for shaded leaf
% --- Set up tridiagonal matrix

%     %---------------------------------------------------------------------
%
%     associate ( &
%                                                % *** Input ***
%     ncan      => mlcanopy_inst%ncan       , &  % Number of aboveground layers
%     nbot      => mlcanopy_inst%nbot       , &  % Index for bottom leaf layer
%     ntop      => mlcanopy_inst%ntop       , &  % Index for top leaf layer
%     dpai      => mlcanopy_inst%dpai       , &  % Layer plant area index (m2/m2)
%     swskyb    => mlcanopy_inst%swskyb     , &  % Atmospheric direct beam solar radiation (W/m2)
%     swskyd    => mlcanopy_inst%swskyd     , &  % Atmospheric diffuse solar radiation (W/m2)
%     fracsun   => mlcanopy_inst%fracsun    , &  % Sunlit fraction of canopy layer
%     fracsha   => mlcanopy_inst%fracsha    , &  % Shaded fraction of canopy layer
%     albsoib   => surfalb_inst%albgrd_col  , &  % Direct beam albedo of ground (soil)
%     albsoid   => surfalb_inst%albgri_col  , &  % Diffuse albedo of ground (soil)
%                                                % *** Output ***
%     swleaf    => mlcanopy_inst%swleaf     , &  % Leaf absorbed solar radiation (W/m2 leaf)
%     swveg     => mlcanopy_inst%swveg      , &  % Absorbed solar radiation, vegetation (W/m2)
%     swvegsun  => mlcanopy_inst%swvegsun   , &  % Absorbed solar radiation, sunlit canopy (W/m2)
%     swvegsha  => mlcanopy_inst%swvegsha   , &  % Absorbed solar radiation, shaded canopy (W/m2)
%     albcan    => mlcanopy_inst%albcan     , &  % Albedo above canopy
%     swsoi     => mlcanopy_inst%swsoi        &  % Absorbed solar radiation, ground (W/m2)
%     )

%---------------------------------------------------------------------
% Set up and solve tridiagonal system of equations for radiative fluxes
%---------------------------------------------------------------------

% --- Set up tridiagonal matrix

for ib = 1:vegetation.params.numrad   % Process each waveband
    for f = 1:vegetation.canopy.num_exposedvegp
        p = vegetation.canopy.filter_exposedvegp(f);
        g = p;
        
        % Zero out radiative fluxes for all layers
        
        swup(p,1,ib) = 0.; %0
        swdn(p,1,ib) = 0.;  %0
        
        for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p)  % ic = 1:ncan
            swup(p,ic,ib) = 0.;
            swdn(p,ic,ib) = 0.;
            vegetation.flux.swleaf(p,ic,isun,ib) = 0.;
            vegetation.flux.swleaf(p,ic,isha,ib) = 0.;
        end
        % There are two equations for each canopy layer and the soil. The first
        % equation is the upward flux and the second equation is the downward flux.
        
        m = 0; % Initialize equation index for tridiagonal matrix
              
%         vegetation.flux.albsoib(p,ib)  = 0.1; %0.3 %surfalb_inst%albgrd_col  ; % Direct beam albedo of ground (soil)
%         vegetation.flux.albsoid(p,ib)  = 0.2; %0.5 %surfalb_inst%albgri_col  ; % Diffuse albedo of ground (soil)

        % Soil: upward flux

        m = m + 1;
        a(m) = 0.;
        b(m) = 1.;
        c(m) = -vegetation.flux.albsoid(g,ib);
        d(m) = vegetation.atmos.swskyb(p,ib) * tbj(p,1) * vegetation.flux.albsoib(g,ib); %tbcum(p,0)  vegetation.flux.albsoib(c,ib)
        
        % Soil: downward flux
        
        refld = (1. - td(p,vegetation.canopy.nbot(p))) * rho(p,ib);
        trand = (1. - td(p,vegetation.canopy.nbot(p))) * tau(p,ib) + td(p,vegetation.canopy.nbot(p));
        aic = refld - trand * trand / refld;
        bic = trand / refld;
        
        m = m + 1;
        a(m) = -aic;
        b(m) = 1.;
        c(m) = -bic;
        d(m) = vegetation.atmos.swskyb(p,ib) * tbj(p,vegetation.canopy.nbot(p)) * (1. - tb(p,vegetation.canopy.nbot(p))) * (tau(p,ib) - rho(p,ib) * bic);
        
        % Leaf layers, excluding top layer
        
        for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)-1
            % Upward flux
            
            refld = (1. - td(p,ic)) * rho(p,ib);
            trand = (1. - td(p,ic)) * tau(p,ib) + td(p,ic);
            fic = refld - trand * trand / refld;
            eic = trand / refld;
            
            m = m + 1;
            a(m) = -eic;
            b(m) = 1.;
            c(m) = -fic;
            d(m) = vegetation.atmos.swskyb(p,ib) * tbj(p,ic) * (1. - tb(p,ic)) * (rho(p,ib) - tau(p,ib) * eic);
            
            % Downward flux
            
            refld = (1. - td(p,ic+1)) * rho(p,ib);
            trand = (1. - td(p,ic+1)) * tau(p,ib) + td(p,ic+1);
            aic = refld - trand * trand / refld;
            bic = trand / refld;
            
            m = m + 1;
            a(m) = -aic;
            b(m) = 1.;
            c(m) = -bic;
            d(m) = vegetation.atmos.swskyb(p,ib) * tbj(p,ic+1) * (1. - tb(p,ic+1)) * (tau(p,ib) - rho(p,ib) * bic);
            
        end
        
        % Top canopy layer: upward flux
        
        ic = vegetation.canopy.ntop(p);
        refld = (1. - td(p,ic)) * rho(p,ib);
        trand = (1. - td(p,ic)) * tau(p,ib) + td(p,ic);
        fic = refld - trand * trand / refld;
        eic = trand / refld;
        
        m = m + 1;
        a(m) = -eic;
        b(m) = 1.;
        c(m) = -fic;
        d(m) = vegetation.atmos.swskyb(p,ib) * tbj(p,ic) * (1 - tb(p,ic)) * (rho(p,ib) - tau(p,ib) * eic); %tau-rho in matlab exercises?
        
        % Top canopy layer: downward flux
        
        m = m + 1;
        a(m) = 0.;
        b(m) = 1.;
        c(m) = 0.;
        d(m) = vegetation.atmos.swskyd(p,ib);
        
        % --- Solve tridiagonal equations for fluxes
        
        [u] = tridiagonal_solver (a, b, c, d, m);
        
        
        %tridiag (atri, btri, ctri, dtri, utri, m)
        
        % Now copy the solution (u) to the upward (swup) and downward (swdn) fluxes for each layer
        % swup - Upward diffuse solar flux above layer
        % swdn - Downward diffuse solar flux onto layer
        
        m = 0;
        
        % Soil fluxes
        
        m = m + 1;
        swup(p,1,ib) = u(m); %0
        m = m + 1;
        swdn(p,1,ib) = u(m); %0
        
        % Leaf layer fluxes
        
        for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
            m = m + 1;
            swup(p,ic,ib) = u(m);
            m = m + 1;
            swdn(p,ic,ib) = u(m);
        end
    end
end

%---------------------------------------------------------------------
% Compute fluxes
%---------------------------------------------------------------------

for ib = 1:vegetation.params.numrad
    for f = 1:vegetation.canopy.num_exposedvegp
        p = vegetation.canopy.filter_exposedvegp(f);
        c = p;
        
        % Solar radiation absorbed by ground (soil)

        swbeam = tbj(p,1) * vegetation.atmos.swskyb(p,ib);  %tbj(p,0)
        swabsb = swbeam * (1. - vegetation.flux.albsoib(c,ib));
        swabsd = swdn(p,1,ib) * (1. - vegetation.flux.albsoid(c,ib));  %swdn(p,0,ib)
        
        vegetation.flux.swdn(p,ib) = swdn(p,1,ib);
        vegetation.flux.swsoi(p,ib) = swabsb + swabsd;
        
        % Leaf layer fluxes
        
        vegetation.flux.swveg(p,ib) = 0.;
        vegetation.flux.swvegsun(p,ib) = 0.;
        vegetation.flux.swvegsha(p,ib) = 0.;
        
        for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
            
            % Downward direct beam incident on layer and absorbed direct
            % beam and diffuse for layer. Note special case for first
            % leaf layer, where the upward flux from below is from the ground.
            % The ground is ic=0, but nbot-1 will not equal 0 if there are lower
            % canopy layers without leaves.
            
            swbeam = tbj(p,ic) * vegetation.atmos.swskyb(p,ib);
            swabsb = swbeam * (1. - tb(p,ic)) * (1. - omega(p,ib));
            if (ic == vegetation.canopy.nbot(p))
                icm1 = 1; %=0
            else 
                icm1 = ic - 1;                    
            end
            
            swabsd = (swdn(p,ic,ib) + swup(p,icm1,ib)) * (1. - td(p,ic)) * (1. - omega(p,ib));
            
            % Absorbed radiation for shaded and sunlit portions of layer
            
            swsha = swabsd * vegetation.flux.fracsha(p,ic);
            swsun = swabsd * vegetation.flux.fracsun(p,ic) + swabsb;
            
            % Per unit sunlit and shaded leaf area
            
            vegetation.flux.swleaf(p,ic,isun,ib) = swsun / (vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic));
            vegetation.flux.swleaf(p,ic,isha,ib) = swsha / (vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic));
            
            % Sum solar radiation absorbed by vegetation and sunlit/shaded leaves
            
            vegetation.flux.swveg(p,ib) = vegetation.flux.swveg(p,ib) + (swabsb + swabsd);
            vegetation.flux.swvegsun(p,ib) = vegetation.flux.swvegsun(p,ib) + swsun;
            vegetation.flux.swvegsha(p,ib) = vegetation.flux.swvegsha(p,ib) + swsha;
            
        end
                
        % Albedo
        
        suminc = vegetation.atmos.swskyb(p,ib) + vegetation.atmos.swskyd(p,ib);
        if (suminc > 0.)
            vegetation.flux.albcan(p,ib) = swup(p,vegetation.canopy.ntop(p),ib) / suminc; % swup(p, vegetation.canopy.ntop(p), ib)
        else 
            vegetation.flux.albcan(p,ib) = 0.;
        end
        
        % Conservation check for total radiation balance: absorbed = incoming - outgoing
        
        sumref = vegetation.flux.albcan(p,ib) * (vegetation.atmos.swskyb(p,ib) + vegetation.atmos.swskyd(p,ib));
        sumabs = suminc - sumref;
        error = sumabs - (vegetation.flux.swveg(p,ib) + vegetation.flux.swsoi(p,ib));
        if (abs(error) > 1.e-03)
            disp ('ERROR: NormanRadiation: total solar conservation error');
        end
        
        % Sunlit and shaded absorption
        
        error = (vegetation.flux.swvegsun(p,ib) + vegetation.flux.swvegsha(p,ib)) - vegetation.flux.swveg(p,ib);
        if (abs(error) > 1.e-03)
            disp ('ERROR: NormanRadiation: sunlit/shade solar conservation error');
        end
        
    end            % end patch loop
end             % end waveband loop

end
