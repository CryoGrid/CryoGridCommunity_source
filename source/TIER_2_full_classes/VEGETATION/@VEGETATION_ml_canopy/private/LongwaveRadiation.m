
function [vegetation] = LongwaveRadiation(vegetation)

% emleaf   = vegetation.params.emleaf            ; %Leaf emissivity (-)
% vegetation.mlcanopyinst.tg       = vegetation.mlcanopyinst.tg         ; %Soil surface temperature (K)
% vegetation.mlcanopyinst.ncan     = vegetation.mlcanopyinst.ncan       ; %Number of aboveground layers
% vegetation.canopy.nbot     = vegetation.canopy.nbot       ; %Index for bottom leaf layer
% vegetation.canopy.ntop     = vegetation.canopy.ntop       ; %Index for top leaf layer
% vegetation.canopy.dpai     = vegetation.canopy.dpai       ; %Layer plant area index (m2/m2)
% vegetation.mlcanopyinst.irsky    = vegetation.mlcanopyinst.irsky      ; %Atmospheric longwave radiation (W/m2)
% vegetation.flux.fracsun  = vegetation.flux.fracsun    ; %Sunlit fraction of canopy layer
% vegetation.flux.fracsha  = vegetation.flux.fracsha    ; %Shaded fraction of canopy layer
% vegetation.mlcanopyinst.tleaf    = vegetation.mlcanopyinst.tleaf      ; %Leaf temperature (K)
% vegetation.flux.td       = vegetation.flux.td         ; %Exponential transmittance of diffuse radiation through a single leaf layer

% *** Output ***
% vegetation.mlcanopyinst.irleaf   = vegetation.mlcanopyinst.irleaf     ; %Leaf absorbed longwave radiation for canopy layer (W/m2 leaf)
% vegetation.flux.irveg    = vegetation.flux.irveg      ; %Absorbed longwave radiation, vegetation (W/m2)
% vegetation.mlcanopyinst.ircan    = vegetation.mlcanopyinst.ircan      ; %Upward longwave radiation above canopy (W/m2)
% vegetation.flux.irsoi    = vegetation.flux.irsoi        ; %Absorbed longwave radiation, ground (W/m2)

% tveg = 273.15 + 25;             % Canopy temperature (K)
% tgrnd = 273.15 + 20;            % Ground temperature (K)

% % % --- paramsmeters
% params.sigma = 5.67e-08;             % Stefan-Boltzmann constant (W/m2/K4)
% params.tfrz = 273.15;                % Freezing point of water (K)
% params.emgrnd = 1.00;                % Ground (soil) emissivity
% emgrnd = 0.96;              % Ground (soil) emissivity

isun = vegetation.params.sun;
isha = vegetation.params.sha;
num_exposedvegp = vegetation.canopy.num_exposedvegp;
filter_exposedvegp = vegetation.canopy.filter_exposedvegp;

for f = 1:num_exposedvegp
    p = filter_exposedvegp(f);
    
    % Zero out radiative fluxes for all layers
    
    irup(p,1) = 0.; %p,0
    irdn(p,1) = 0.; %p,0
    
    for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran 1:ncan
        irup(p,ic) = 0.0;
        irdn(p,ic) = 0.0;
        vegetation.mlcanopyinst.irleaf(p,ic) = 0.0;
    end
    
    % Ground (soil) emissivity
    
    %CHANGE SEBASTIAN
    %emg = 0.96; %SET AVERAGE EMISSIVITY of the GROUND 
    emg = vegetation.flux.emissivity_ground;
    
    %------------------------------------------------------------------
    % Leaf scattering coefficient and terms for longwave radiation reflected
    % and transmitted by a layer
    %------------------------------------------------------------------
    
    omega = 1. - vegetation.pftcon.emleaf(p);
    
    % Intercepted radiation is reflected
    
    rho = omega ;
    tau = 0.;
%     
    % Intercepted radiation is both reflected and transmitted
    
%          rho = omega * 0.5;
%          tau = omega * 0.5;
    
    %------------------------------------------------------------------
    % Emitted longwave radiation is weighted average of sunlit and shaded leaves
    %------------------------------------------------------------------
    sb = vegetation.physcon.sigma; % Stefan boltzman const
    for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
%         vegetation.flux.td(p,ic) = 0.1; %where should this be created?!! Exponetial transmittance
        ir_source_sun = vegetation.pftcon.emleaf(p) * sb * vegetation.mlcanopyinst.tleaf(p,ic,isun)^4;
        ir_source_sha = vegetation.pftcon.emleaf(p) * sb * vegetation.mlcanopyinst.tleaf(p,ic,isha)^4;
        ir_source(ic) = (ir_source_sun * vegetation.flux.fracsun(p,ic) + ir_source_sha * vegetation.flux.fracsha(p,ic)) * (1. - vegetation.flux.td(p,ic));
    end
    
    %------------------------------------------------------------------
    % Set up and solve tridiagonal system of equations for upward and forwnward fluxes
    %------------------------------------------------------------------
    
    % There are two equations for each leaf layer and the soil. The first
    % equation is the upward flux and the second equation is the downward flux.
    
    m = 0;
    
    % Soil: upward flux
    
    m = m + 1;
    atri(m) = 0.;
    btri(m) = 1.;
    ctri(m) = -(1. - emg);
    %dtri(m) = emg * sb * vegetation.mlcanopyinst.tg(p)^4;  %SET outgoing radiation here, but only emitted part!
    dtri(m) = vegetation.flux.Lout_ground;
    
    % Soil: downward flux
    
    refld = (1. - vegetation.flux.td(p,vegetation.canopy.nbot(p))) * rho;
    trand = (1. - vegetation.flux.td(p,vegetation.canopy.nbot(p))) * tau + vegetation.flux.td(p,vegetation.canopy.nbot(p));
    aic = refld - trand * trand / refld;
    bic = trand / refld;
    
    m = m + 1;
    atri(m) = -aic;
    btri(m) = 1.;
    ctri(m) = -bic;
    dtri(m) = (1. - bic) * ir_source(vegetation.canopy.nbot(p));
    
    % Leaf layers, excluding top layer
    
    for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)-1
        
        % Upward flux
        
        refld = (1. - vegetation.flux.td(p,ic)) * rho;
        trand = (1. - vegetation.flux.td(p,ic)) * tau + vegetation.flux.td(p,ic);
        fic = refld - trand * trand / refld;
        eic = trand / refld;
        
        m = m + 1;
        atri(m) = -eic;
        btri(m) = 1.;
        ctri(m) = -fic;
        dtri(m) = (1. - eic) * ir_source(ic);
        
        % Downward flux
        
        refld = (1. - vegetation.flux.td(p,ic+1)) * rho;
        trand = (1. - vegetation.flux.td(p,ic+1)) * tau + vegetation.flux.td(p,ic+1);
        aic = refld - trand * trand / refld;
        bic = trand / refld;
        
        m = m + 1;
        atri(m) = -aic;
        btri(m) = 1.;
        ctri(m) = -bic;
        dtri(m) = (1. - bic) * ir_source(ic+1);
        
    end
    
    % Top canopy layer: upward flux
    
    ic = vegetation.canopy.ntop(p);
    refld = (1. - vegetation.flux.td(p,ic)) * rho;
    trand = (1. - vegetation.flux.td(p,ic)) * tau + vegetation.flux.td(p,ic);
    fic = refld - trand * trand / refld;
    eic = trand / refld;
    
    m = m + 1;
    atri(m) = -eic;
    btri(m) = 1.;
    ctri(m) = -fic;
    dtri(m) = (1. - eic) * ir_source(ic);
    
    % Top canopy layer: downward flux
    
    m = m + 1;
    atri(m) = 0.;
    btri(m) = 1.;
    ctri(m) = 0.;
    dtri(m) = vegetation.mlcanopyinst.irsky(p);
    
    % Solve tridiagonal system of equations for upward and downward fluxes
    
    [utri] = tridiagonal_solver (atri, btri, ctri, dtri, m);
    
    % Now copy the solution (utri) to the upward (irup) and downward (irdn)
    % fluxes for each layer
    % irup =  Upward longwave flux above layer
    % irdn =  Downward longwave flux onto layer
    
    m = 0;
    
    % Soil fluxes
    
    m = m + 1;
    irup(p,1) = utri(m); %(p,0)
    m = m + 1;
    irdn(p,1) = utri(m); %(p,0)
    
    % Leaf layer fluxes
    
    for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
        m = m + 1;
        irup(p,ic) = utri(m);
        m = m + 1;
        irdn(p,ic) = utri(m);
    end
    
    %------------------------------------------------------------------
    % Compute fluxes
    %------------------------------------------------------------------
    
    % Absorbed longwave radiation for ground (soil)
    
    vegetation.flux.irsoi(p) = irdn(p,1) - irup(p,1); %(p,0)%(p,0)
    vegetation.flux.irdn = irdn(p,1);
    % Leaf layer fluxes
    
    vegetation.flux.irveg(p) = 0.;
    
    for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
        
        % Absorbed longwave radiation for layer. Note special case for first
        % leaf layer, where the upward flux from below is from the ground.
        % The ground is ic=0, but vegetation.canopy.nbot-1 will not equal 0 if there are lower
        % canopy layers without leaves.
        
        if (ic == vegetation.canopy.nbot(p))
            icm1 = 1; %0
        else
            icm1 = ic-1; %ic-1
        end
        irabs = vegetation.pftcon.emleaf(p) * (irdn(p,ic)+irup(p,icm1)) * (1. - vegetation.flux.td(p,ic)) - 2. * ir_source(ic);
        vegetation.mlcanopyinst.irleaf(p,ic) = irabs / vegetation.canopy.dpai(p,ic);
        
        % Sum longwave radiation absorbed by vegetation
        
        vegetation.flux.irveg(p) = vegetation.flux.irveg(p) + irabs;
        
    end
    
    % Canopy emitted longwave radiation
    
    vegetation.mlcanopyinst.ircan(p) = irup(p,vegetation.canopy.ntop(p));
    
    %------------------------------------------------------------------
    % Conservation check
    %------------------------------------------------------------------
    
    % Total radiation balance: absorbed = incoming - outgoing
    
    sumabs = vegetation.mlcanopyinst.irsky(p) - vegetation.mlcanopyinst.ircan(p);
    error = sumabs - (vegetation.flux.irveg(p) + vegetation.flux.irsoi(p));
    if (abs(error) > 1.e-03)
        disp('ERROR: LongwaveRadiationMod: total longwave conservation error');
    end
end
end



%-----------------------------------------------------------------------
% %DESCRIPTION:
% Calculate longwave radiation transfer through canopy

% %DESCRIPTION:
% Longwave radiation transfer through canopy using Norman (1979)
%
% %USES:


%     use decompMod, only : bounds_type
%     use pftconMod, only : pftcon
%     use PatchType, only : patch
%     use MathToolsMod, only : tridiag
%     use CanopyFluxesMultilayerType, only : mlcanopy_type

%     % %ARGUMENTS:
%     implicit none
%     type(bounds_type), intent(in) :: bounds
%     integer, intent(in) :: num_exposedvegp        % Number of non-snow-covered veg points in CLM patch filter
%     integer, intent(in) :: filter_exposedvegp(:)  % CLM patch filter for non-snow-covered vegetation
%     type(mlcanopy_type), intent(inout) :: mlcanopy_inst
%     %
%     % %LOCAL VARIABLES:
%     integer  :: f                             % Filter index
%     integer  :: p                             % Patch index for CLM g/l/c/p hierarchy
%     integer  :: ic                            % Aboveground layer index
%     integer  :: icm1                          % Layer below ic (ic-1)
%     real(r8) :: emg                           % Ground (soil) emissivity
%     real(r8) :: sumabs                        % Absorbed radiation for energy conservation check
%     real(r8) :: error                         % Error check
%     real(r8) :: omega                         % Leaf scattering coefficient
%     real(r8) :: rho                           % Leaf reflectance
%     real(r8) :: tau                           % Leaf transmittance
%     real(r8) :: trand                         % Term for longwave radiation transmitted by layer
%     real(r8) :: refld                         % Term for longwave radiation reflected by layer
%     real(r8) :: ir_source_sun                 % Longwave radiation emitted by sunlit leaf (W/m2)
%     real(r8) :: ir_source_sha                 % Longwave radiation emitted by shaded leaf (W/m2)
%     real(r8) :: ir_source(nlevcan)            % Longwave radiation emitted by leaf layer (W/m2)
%     integer  :: m                             % Index to the tridiagonal matrix
%     real(r8) :: aic, bic                      % Intermediate terms for tridiagonal matrix
%     real(r8) :: eic, fic                      % Intermediate terms for tridiagonal matrix
%     integer, parameter :: neq = (nlevcan+1)*2 % Number of tridiagonal equations to solve
%     real(r8) :: atri(neq), btri(neq)          % Entries in tridiagonal matrix
%     real(r8) :: ctri(neq), dtri(neq)          % Entries in tridiagonal matrix
%     real(r8) :: utri(neq)                     % Tridiagonal solution
%     real(r8) :: irabs                         % Absorbed longwave flux (W/m2 ground)
%     real(r8) :: irup(bounds%begp:bounds%endp,0:nlevcan) % Upward longwave flux above canopy layer (W/m2 ground)
%     real(r8) :: irdn(bounds%begp:bounds%endp,0:nlevcan) % Downward longwave flux onto canopy layer (W/m2 ground)
%---------------------------------------------------------------------


    