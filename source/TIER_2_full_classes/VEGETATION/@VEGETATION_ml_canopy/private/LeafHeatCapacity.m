%-----------------------------------------------------------------------
function [vegetation] = LeafHeatCapacity(vegetation, p, ic) 
            
%(p, ic, mlcanopy_inst, clm_varcon, clm_varctl, pftconMod, PatchType, CanopyFluxesMultilayerType, pftcon)
%
% %DESCRIPTION:
% Leaf heat capacity
%
% %USES:
% cpliq = vegetation.pftcon.cpliq; %use clm_varcon
% cpbio = vegetation.pftcon.cpbio; %use clm_varcon
% no_storage = vegetation.pftcon.no_storage; %use clm_varctl
% dpai = vegetation.canopy.dpai;
% p = vegetation.mlcanopyinst.p;            % Patch index for CLM g/l/c/p hierarchy
% ic = 1;%vegetation.mlcanopyinst.ic;                % Aboveground layer index

% mlcanopy_inst

% %LOCAL VARIABLES:
% lma                          % Leaf carbon mass per area (kg C / m2 leaf)
% dry_weight                   % Leaf dry mass per area (kg DM / m2 leaf)
% fresh_weight                 % Leaf fresh mass per area (kg FM / m2 leaf)
% leaf_water                   % Leaf water (kg H2O / m2 leaf)
% fcarbon = vegetation.pftcon.fcarbon;  % Fraction of dry biomass that is carbon
% fwater = vegetation.pftcon.fwater;   % Fraction of fresh biomass that is water
%---------------------------------------------------------------------

% *** Input ***
% slatop = vegetation.pftcon.slatop;           
% slatop = pftcon.slatop           &  % Specific leaf area at top of canopy (m2/gC)
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Leaf heat capacity - need to convert specific leaf area (m2/gC) to
% leaf mass per area (kgC/m2) and then convert to dry weight (assume
% carbon is 50% of dry biomass). Then need to convert dry biomass to
% fresh biomass (assume 70% of fresh biomass is water). Then remember
% that 70% of fresh biomass is water when calculating heat capacity.

% vegetation.mlcanopyinst.cpleaf = ones(1,vegetation.mlcanopyinst.ncan); 

if (vegetation.canopy.dpai(p,ic) > 0) % >0 leaf layer
    lma = 1 ./ vegetation.pftcon.slatop(p) .* 0.001;                                % m2 / g C -> kg C / m2
    dry_weight = lma ./ vegetation.pftcon.fcarbon;                                 % kg C / m2 -> kg DM / m2
    fresh_weight = dry_weight ./ (1. - vegetation.pftcon.fwater);                      % kg DM / m2 -> kg FM / m2
    leaf_water = vegetation.pftcon.fwater .* fresh_weight;                        % Leaf water (kg H2O / m2 leaf)
    vegetation.mlcanopyinst.cpleaf(p,ic) = vegetation.pftcon.cpbio .* dry_weight + vegetation.pftcon.cpliq .* leaf_water;        % Heat capacity (J/K/m2 leaf)
else
    vegetation.mlcanopyinst.cpleaf(p,ic) = 0.;
end

    
%     % Use very low leaf heat capacity to reduce storage term
%     no_storage = 1; 
%     if no_storage == 1 
%         vegetation.mlcanopyinst.cpleaf(p,ic) = 3350; %bonan book value
%     end
%     
    % Output
%     vegetation.mlcanopyinst.cpleaf = cpleaf;
    
end