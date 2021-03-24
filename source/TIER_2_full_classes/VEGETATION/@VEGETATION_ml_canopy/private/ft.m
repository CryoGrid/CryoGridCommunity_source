%-----------------------------------------------------------------------
function [pstr] = ft (vegetation, t, ha)

% DESCRIPTION:
% Photosynthesis temperature response

% USES:
% rgasc = vegetation.physcon.rgasc;
% tfrz = vegetation.physcon.tfrz;
% t = vegetation.mlcanopyinst.tleaf(p,ic,il); 
% ha = vegetation.leaf.kcha; 

 
% ha = vegetation.mlcanopyinst.ha; %dummy values for now

%     ARGUMENTS:
%     implicit none
%     real(r8) :: tl                   % Leaf temperature (K)
%     real(r8) :: ha                   % Activation energy (J/mol)
%
%     LOCAL VARIABLES:
%     real(r8) :: ans                  % Temperature function value
%---------------------------------------------------------------------

pstr = exp(ha / (vegetation.physcon.rgasc * (vegetation.physcon.tfrz+25)) * (1 - (vegetation.physcon.tfrz+25) / t) );

end