
%-----------------------------------------------------------------------
function [ftv] = fth (vegetation, t, hd, se, c)
%(tl, hd, se, c, clm_varcon)

% DESCRIPTION:
% Photosynthesis temperature inhibition
% tl = vegetation.mlcanopyinst.t1; %dummy values for now
% hd = vegetation.mlcanopyinst.hd; %dummy values for now
% se = vegetation.mlcanopyinst.se; %dummy values for now
% c = vegetation.mlcanopyinst.c; %dummy values for now

% ARGUMENTS:
%     implicit none
%     real(r8) :: tl                   % Leaf temperature (K)
%     real(r8) :: hd                   % Deactivation energy (J/mol)
%     real(r8) :: se                   % Entropy term (J/mol/K)
%     real(r8) :: c                    % Scaling factor for high temperature inhibition (25 C = 1.0)
%
% %LOCAL VARIABLES:
%     real(r8) :: ans                  % Temperature function value
%---------------------------------------------------------------------

ftv = c ./ (1+ exp((-hd + se*t) ./ (vegetation.physcon.rgasc*t)));


end