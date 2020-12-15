% module CanopyNitrogenProfileMod

%-----------------------------------------------------------------------
% DESCRIPTION:
% Canopy profile of nitrogen and photosynthetic capacity

% USES:
%   shr_kind_r8 ; shr_kind_mod

% PUBLIC TYPES:

% %PUBLIC MEMBER FUNCTIONS:
%   public :: CanopyNitrogenProfile
%-----------------------------------------------------------------------

function [vegetation] = CanopyNitrogenProfile(vegetation) %, p, ic
%(num_exposedvegp, filter_exposedvegp, mlcanopy_inst, vcmaxpft, c3psn, clm_varct1, tfrz, patch, pftconMod)

% DESCRIPTION:
% Canopy profile of nitrogen and photosynthetic capacity

% USES:
% use_clm45kn = vegetation.physcon.use_clm45kn; % clm_varctl
% use_acclim = vegetation.physcon.use_acclim; % clm_varctl
tfrz = vegetation.physcon.tfrz;                          % Freezing point of water (K)

% ARGUMENTS:
num_exposedvegp = vegetation.canopy.num_exposedvegp;       % Number of non-snow-covered veg points in CLM patch filter
filter_exposedvegp = vegetation.canopy.filter_exposedvegp;  % CLM patch filter for non-snow-covered vegetation

% %LOCAL VARIABLES:
%     f               % Filter index
%     p               % Patch index for CLM g/l/c/p hierarchy
%     ic              % Aboveground layer index
%     vcmax25top      % Canopy top - Maximum carboxylation rate at 25C (umol/m2/s)
%     jmax25top       % Canopy top - C3: Maximum electron transport rate at 25C (umol/m2/s)
%     rd25top         % Canopy top - Leaf respiration rate at 25C (umol CO2/m2/s)
%     kp25top         % Canopy top - C4: Initial slope of CO2 response curve at 25C (mol/m2/s)
%     kn              % Leaf nitrogen decay coefficient
%     nscale          % Nitrogen scaling coefficient
%---------------------------------------------------------------------
% Input
% vcmaxpft    = pftcon.vcmaxpft;           % Maximum carboxylation rate at 25C (umol/m2/s)
% c3psn       = pftcon.c3psn    ;          % Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant
% vcmaxpft    = vegetation.pftcon.vcmaxpft;
% c3psn       = vegetation.pftcon.c3psn;

% tacclim     = vegetation.mlcanopyinst.tacclim  ;       % Average air temperature for acclimation (K)
% ncan        = vegetation.mlcanopyinst.ncan  ;      % Number of aboveground layers
% dpai        = vegetation.canopy.dpai   ;     % Layer plant area index (m2/m2)
% sumpai      = vegetation.mlcanopyinst.sumpai  ;    % Cumulative plant area index (m2/m2)

% Output
% vcmax25     = vegetation.leaf.vcmax25  ;   % Leaf maximum carboxylation rate at 25C for canopy layer (umol/m2/s)
% jmax25      = vegetation.leaf.jmax25    ;  % C3 - maximum electron transport rate at 25C for canopy layer (umol/m2/s)
% rd25        = vegetation.leaf.rd25       ; % Leaf respiration rate at 25C for canopy layer (umol CO2/m2/s)
% kp25        = vegetation.leaf.kp25_c4    ;% C4 - initial slope of CO2 response curve at 25C for canopy layer (mol/m2/s)

%     for f = 1:num_exposedvegp
%        p = filter_exposedvegp(f);

for f = 1:num_exposedvegp
    p = filter_exposedvegp(f);
    % vcmax and other parameters (at 25C and top of canopy). jmax acclimation
    % from Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190
    
    vcmax25top = vegetation.pftcon.vcmaxpft(p);
    
    % % %     if (round(c3psn(p)) == 1)
    %         if (use_acclim)
    jmax25top = (2.59  - 0.035 .* min(max((vegetation.mlcanopyinst.tacclim(p)-tfrz),11.),35.)) .* vcmax25top;
    %         else
    %             jmax25top = 1.67 .* vcmax25top;
    %         end
    rd25top = 0.015 .* vcmax25top;
    kp25top = 0.;
    % % %     else
    % % %         jmax25top = 0;
    % % %         rd25top = 0.025 .* vcmax25top;
    % % %         kp25top = 0.02 .* vcmax25top;
    % % %     end
    
    % Leaf nitrogen decay coefficient
    
    %     if (use_clm45kn)
    %         kn = 0.3;
    %     else
    kn = exp(0.00963  .* vcmax25top - 2.43);
    %     end
    
    % Layer values
    
    for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
        if (vegetation.canopy.dpai(p,ic) > 0)   % leaf layer
            nscale = exp (-kn .* vegetation.mlcanopyinst.sumpai(p,ic));
            vegetation.leaf.vcmax25(p,ic) = vcmax25top .* nscale;
            vegetation.leaf.jmax25(p,ic) = jmax25top .* nscale;
            vegetation.leaf.rd25(p,ic) = rd25top .* nscale;
            vegetation.leaf.kp25(p,ic) = kp25top .* nscale;
        else % non-leaf leaf layer
            vegetation.leaf.vcmax25(p,ic) = 0;
            vegetation.leaf.jmax25(p,ic) = 0;
            vegetation.leaf.rd25(p,ic) = 0;
            vegetation.leaf.kp25(p,ic) = 0;
        end
    end
end

% vegetation.leaf.vcmax25     = vcmax25;   % Leaf maximum carboxylation rate at 25C for canopy layer (umol/m2/s)
% vegetation.leaf.jmax25      = jmax25;  % C3 - maximum electron transport rate at 25C for canopy layer (umol/m2/s)
% vegetation.leaf.rd25        = rd25; % Leaf respiration rate at 25C for canopy layer (umol CO2/m2/s)
% vegetation.leaf.kp25        = kp25;% C4 - initial slope of CO2 response curve at 25C for canopy layer (mol/m2/s)


end