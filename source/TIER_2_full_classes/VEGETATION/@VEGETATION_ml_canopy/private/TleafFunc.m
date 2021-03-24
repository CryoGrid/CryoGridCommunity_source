function [vegetation, tleaf_dif] = TleafFunc (vegetation,p,ic,il,tleaf_val)
%(p, ic, il, mlcanopy_inst, tleaf_val, shr_kind, abortutils, clm_varctl,CanopyFluxesMultilayerType, LeafPhotosynthesisMod, LeafTemperatureMod)
% DESCRIPTION:
% Calculate leaf fluxes for an input leaf temperature (tleaf_val) and
% compare the new temperature to the prior temperature. This function
% equals zero when tleaf does not change between iterations.

% USES:
%   iulog = 1; % clm_varctl

% leaf_temp_iter = vegetation.physcon.leaf_temp_iter; % clm_varctl
% p           = vegetation.mlcanopyinst.p;      % Patch index for CLM g/l/c/p hierarchy
% ic          = 1;%vegetation.mlcanopyinst.ic;     % Canopy layer index
% il          = 1;%vegetation.mlcanopyinst.il;     % Sunlit (1) or shaded (2) leaf index

% ARGUMENTS:
% tleaf = vegetation.mlcanopyinst.tleaf;  % Input value for tleaf

% LOCAL VARIABLES:
%   tleaf_dif = 1;              % Difference in tleaf
%---------------------------------------------------------------------

%   associate (tleaf => mlcanopy_inst.tleaf)  % Leaf temperature (K)

if (tleaf_val < 0)
    disp ('TleafFunc error')
end

vegetation.mlcanopyinst.tleaf(p,ic,il) = tleaf_val;

% Leaf photosynthesis and stomatal conductance

[vegetation] = LeafPhotosynthesis(vegetation,p,ic,il);

% Leaf temperature and energy fluxes

[vegetation] = LeafTemperatureMod(vegetation,p,ic,il);

% --- Compare with prior value for leaf temperature

tleaf_dif = vegetation.mlcanopyinst.tleaf(p,ic,il) - tleaf_val;

end