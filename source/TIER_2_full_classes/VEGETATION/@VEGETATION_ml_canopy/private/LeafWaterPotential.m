function [vegetation, leafwp] = LeafWaterPotential (vegetation,p,ic,il,leafwp)
% %ARGUMENTS:
% p              % Patch index for CLM g/l/c/p hierarchy
% ic             % Aboveground layer index
% il             % Sunlit (1) or shaded (2) leaf index
% leafwp         % Leaf water potential (MPa)
% mlcanopy_inst

% %LOCAL VARIABLES:
% dtime                         % Model time step (s)
% y0                            % Leaf water potential at beginning of timestep (MPa)
% dy                            % Change in leaf water potential (MPa)
% a, b                          % Intermediate calculation

vegetation.physcon.head = vegetation.physcon.denh20*vegetation.physcon.grav*1.e-06;
% head                = vegetation.physcon.head;
% capac               = vegetation.pftcon.capac(p);                  % Plant capacitance (mmol H2O/m2 leaf area/MPa)
% psis                = vegetation.mlcanopyinst.psis(p);          % Weighted soil water potential (MPa)
% zs                  = vegetation.mlcanopyinst.zs(p,ic);         % Canopy height for scalar concentration and source (m)
% % dpai                = vegetation.canopy.dpai(p,ic);       % Layer plant area index (m2/m2)
% lsc                 = vegetation.mlcanopyinst.lsc(p,ic);        % Leaf-specific conductance for canopy layer (mmol H2O/m2 leaf/s/MPa)
% trleaf              = vegetation.mlcanopyinst.trleaf(p,ic,il);  % Leaf transpiration flux (mol H2O/m2 leaf/s)
% leafwp              = vegetation.mlcanopyinst.leafwp(p,ic,il);

% Get step size
dtime = vegetation.params.dtime_sub;

% Change in leaf water potential is: dy / dt = (a - y) / b. The integrated change
% over a full model timestep is: dy = (a - y0) * (1 - exp(-dt/b))

if (vegetation.canopy.dpai(p,ic) > 0)  % leaf layer
    y0 = leafwp;
    a = vegetation.mlcanopyinst.psis(p) - vegetation.physcon.head .* vegetation.mlcanopyinst.zs(p,ic) - 1000. .* vegetation.mlcanopyinst.trleaf(p,ic,il) ./ vegetation.mlcanopyinst.lsc(p,ic);
    b = vegetation.leaf.capac(p) ./ vegetation.mlcanopyinst.lsc(p,ic);
    dy = (a-y0) .* (1.-exp(-dtime./b));
    leafwp = y0 + dy;
else % non-leaf layer
    leafwp = 0;
end
%vegetation.mlcanopyinst.leafwp = leafwp;

end


