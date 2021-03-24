%-----------------------------------------------------------------------
function [fth25F] = fth25_g(ground, hd, se)
%
% %DESCRIPTION:
% Scaling factor for photosynthesis temperature inhibition

% hd = ground.STATVAR.vegetation.mlcanopyinst.hd; %dummy values for now
% se = ground.STATVAR.vegetation.mlcanopyinst.se; %dummy values for now

% ARGUMENTS:
%     hd = 1;  % Deactivation energy (J/mol)
%     se = 1;  % Entropy term (J/mol/K)

% LOCAL VARIABLES:
%ans                  % Temperature function value
%---------------------------------------------------------------------

fth25F = 1.  + exp((-hd + se .* (ground.STATVAR.vegetation.physcon.tfrz+25)) ./ (ground.STATVAR.vegetation.physcon.rgasc .* (ground.STATVAR.vegetation.physcon.tfrz+25)));

% ground.STATVAR.vegetation.mlcanopyinst.fth25F = fth25F;
end
