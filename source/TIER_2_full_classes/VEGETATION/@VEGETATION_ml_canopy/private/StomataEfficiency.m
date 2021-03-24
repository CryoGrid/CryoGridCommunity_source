  %-----------------------------------------------------------------------
  function [vegetation, val] = StomataEfficiency (vegetation,p,ic,il,gs_val)
    %, gs_val
    % %DESCRIPTION:
    % Stomata water-use efficiency check and cavitation check to determine maximum gs. 
    % For the stomatal conductance gs_val, calculate photosynthesis and leaf
    % water potential for an increase in stomatal conductance equal to "delta".
    % The returned value is positive if this increase produces a change in
    % photosynthesis > iota*vpd*delta or if the leaf water potential is > minlwp.
    % The returned value is negative if the increase produces a change in
    % photosynthesis < iota*vpd*delta or if the leaf water potential is < minlwp. 
    %
%     % %USES:
%     use pftconMod, only : pftcon
%     use PatchType, only : patch
%     use CanopyFluxesMultilayerType, only : mlcanopy_type
%     %
%     % %ARGUMENTS:
%     implicit none
%     integer, intent(in)  :: p        % Patch index for CLM g/l/c/p hierarchy
% p           = vegetation.mlcanopyinst.p;      % Patch index for CLM g/l/c/p hierarchy
% ic          = vegetation.mlcanopyinst.ic;     % Canopy layer index
% il          = vegetation.mlcanopyinst.il;     % Sunlit (1) or shaded (2) leaf index
% gs_val = 1;  % Value for gs to use in calculations
% gs_val = 1;

%     % %LOCAL VARIABLES:
%     real(r8) :: delta                % Small difference for gs (mol H2O/m2/s)
%     real(r8) :: leafwp               % Current leaf water potential (MPa)
%     real(r8) :: gs2                  % Lower value for gs (mol H2O/m2/s)
%     real(r8) :: an2                  % Leaf photosynthesis at gs2 (umol CO2/m2/s)
%     real(r8) :: gs1                  % Higher value for gs (mol H2O/m2/s)
%     real(r8) :: an1                  % Leaf photosynthesis at gs1 (umol CO2/m2/s)
%     real(r8) :: wue                  % Water-use efficiency check
%     real(r8) :: minpsi               % Cavitation check
%     real(r8) :: val                  % Returned minimum of the two checks
%     %---------------------------------------------------------------------
%     minlwp    = vegetation.pftcon.minlwp ; %pftcon.minlwp      ;  % Minimum leaf water potential (MPa)
%     iota      = vegetation.leaf.iota; %pftcon.iota        ;  % Stomatal water-use efficiency (umol CO2/ mol H2O)
%     psil      = vegetation.mlcanopyinst.psil ;  % Leaf water potential (MPa)
%     an        = vegetation.mlcanopyinst.an   ;  % Leaf net photosynthesis (umol CO2/m2 leaf/s)
%     vpd       = vegetation.mlcanopyinst.vpd  ;  % Leaf vapor pressure deficit (Pa)
%     pref      = vegetation.mlcanopyinst.pref   ;  % Air pressure at reference height (Pa)

    % Specify "delta" as a small difference in gs (mol H2O/m2/s)

    delta = 0.001;
    % Photosynthesis at lower gs (gs_val - delta)

    leafwp = vegetation.mlcanopyinst.psil(p,ic,il);
    gs2 = gs_val - delta;
    [vegetation, leafwp] = StomataFluxes (vegetation, p,ic,il, gs2, leafwp);
    an2 = vegetation.mlcanopyinst.an(p,ic,il);

    % Photosynthesis at higher gs (gs_val)

    leafwp = vegetation.mlcanopyinst.psil(p,ic,il);
    gs1 = gs_val;
    [vegetation, leafwp] = StomataFluxes (vegetation, p,ic,il, gs1, leafwp);
    an1 = vegetation.mlcanopyinst.an(p,ic,il); 

    % Efficiency check: wue < 0 when d(An) / d(gs) < iota * vpd

    wue = (an1 - an2) - vegetation.leaf.iota(p) * delta * (vegetation.mlcanopyinst.vpd(p,ic,il) / vegetation.mlcanopyinst.pref(p));

    % Cavitation check: minpsi < 0 when leafwp < minlwp

    minpsi = leafwp - vegetation.leaf.minlwp(p);

    % Return the minimum of the two checks

    val = min(wue, minpsi);

  end

