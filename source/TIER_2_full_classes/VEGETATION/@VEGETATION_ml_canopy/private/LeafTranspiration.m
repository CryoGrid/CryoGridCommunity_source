function [vegetation] = LeafTranspiration (vegetation,p,ic,il)
%
% %DESCRIPTION:
% Calculate leaf water transpiration flux
% %ARGUMENTS:
%     p            % Patch index for CLM g/l/c/p hierarchy
%     ic           % Aboveground layer index
%     il           % Sunlit (1) or shaded (2) leaf index
%     mlcanopy_inst

% LOCAL VARIABLES:
%     esat                    % Saturation vapor pressure (Pa)
%     desat                   % Temperature derivative of saturation vapor pressure (Pa/K)
%     lambda                  % Latent heat of vaporization (J/mol)
%     gleaf                   % Leaf conductance for transpiration (mol H2O/m2 leaf/s)
%---------------------------------------------------------------------

% p         = vegetation.mlcanopyinst.p;      % Patch index for CLM g/l/c/p hierarchy
% ic        = vegetation.mlcanopyinst.ic;     % Canopy layer index
% il        = vegetation.mlcanopyinst.il;     % Sunlit (1) or shaded (2) leaf index
% dpai      = vegetation.canopy.dpai;   % Layer plant area index (m2/m2)
% tref      = vegetation.mlcanopyinst.tref;   % Air temperature at reference height (K)
% pref      = vegetation.mlcanopyinst.pref;   % Air pressure at reference height (Pa)
% % eair      = vegetation.mlcanopyinst.eair;   % Vapor pressure profile (Pa)
% fdry      = vegetation.mlcanopyinst.fdry;   % Fraction of plant area index that is green and dry
% % tleaf     = vegetation.mlcanopyinst.tleaf;  % Leaf temperature (K)
% gs        = vegetation.mlcanopyinst.gs;     % Leaf stomatal conductance (mol H2O/m2 leaf/s)
% gbv       = vegetation.mlcanopyinst.gbv;    % Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
% % Output
% trleaf    = vegetation.mlcanopyinst.trleaf; % Leaf transpiration flux (mol H2O/m2 leaf/s)

if (vegetation.canopy.dpai(p,ic) > 0.)   % leaf layer
    
    % Saturation vapor pressure
    
    [esat, ~] = Satvap (vegetation.mlcanopyinst.tleaf(p,ic,il));
    
    % Latent heat of vaporization
    
    [lambda] = LatVap (vegetation.mlcanopyinst.tref(p),vegetation);
    
    % Leaf conductance for transpiration
    
    gleaf = vegetation.mlcanopyinst.gs(p,ic,il) .* vegetation.mlcanopyinst.gbv(p,ic,il) ./ (vegetation.mlcanopyinst.gs(p,ic,il) + vegetation.mlcanopyinst.gbv(p,ic,il));
    
    % Transpiration flux: mol H2O/m2/s
    
    vegetation.mlcanopyinst.trleaf(p,ic,il) = gleaf .* vegetation.mlcanopyinst.fdry(p,ic) .* (esat - vegetation.mlcanopyinst.eair(p,ic)) ./ vegetation.mlcanopyinst.pref(p);
    
else
    
    vegetation.mlcanopyinst.trleaf(p,ic,il) = 0.;
    
end

% Output
% vegetation.mlcanopyinst.trleaf    = trleaf      ;  % Leaf transpiration flux (mol H2O/m2 leaf/s)

end
