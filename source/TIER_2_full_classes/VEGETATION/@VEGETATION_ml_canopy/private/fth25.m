%-----------------------------------------------------------------------
function [fth25F] = fth25(vegetation, hd, se)
%
% %DESCRIPTION:
% Scaling factor for photosynthesis temperature inhibition

% hd = vegetation.mlcanopyinst.hd; %dummy values for now
% se = vegetation.mlcanopyinst.se; %dummy values for now

% ARGUMENTS:
%     hd = 1;  % Deactivation energy (J/mol)
%     se = 1;  % Entropy term (J/mol/K)

% LOCAL VARIABLES:
%ans                  % Temperature function value
%---------------------------------------------------------------------

fth25F = 1.  + exp((-hd + se .* (vegetation.physcon.tfrz+25)) ./ (vegetation.physcon.rgasc .* (vegetation.physcon.tfrz+25)));

% vegetation.mlcanopyinst.fth25F = fth25F;
end
