%-----------------------------------------------------------------------
function [vegetation, leafwp] = StomataFluxes (vegetation,p,ic,il,gs_val,leafwp)

% DESCRIPTION:
% Calculate leaf temperature, energy fluxes, photosynthesis, and leaf
% water potential for a specificed stomatal conductance (gs_val)

% leaf_temp_iter = vegetation.physcon.leaf_temp_iter;

%     % ARGUMENTS:
%     leafwp    % Leaf water potential (MPa)
%     gs = vegetation.mlcanopyinst.gs;
%---------------------------------------------------------------------

% gs_val(p,ic,il) = vegetation.mlcanopyinst.gs(p,ic,il);  % Leaf stomatal conductance (mol H2O/m2 leaf/s)

% Use specified gs (gs_val)
vegetation.mlcanopyinst.gs(p,ic,il) = gs_val;

% vegetation.mlcanopyinst.gs = gs;
% Leaf temperature and energy fluxes

%if %(vegetation.physcon.leaf_temp_iter == 1)
[vegetation] = LeafTemperatureMod (vegetation,p,ic,il); 
%end

% Leaf photosynthesis

[vegetation] = LeafPhotosynthesis (vegetation,p,ic,il);

% Leaf transpiration

% if (~ leaf_temp_iter)
[vegetation] = LeafTranspiration (vegetation,p,ic,il);
% end

% Leaf water potential
[vegetation, leafwp] = LeafWaterPotential(vegetation,p,ic,il,leafwp);

end


  