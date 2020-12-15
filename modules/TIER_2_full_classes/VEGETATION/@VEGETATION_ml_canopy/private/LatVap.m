  function [lambda] = LatVap (t, vegetation)
%DESCRIPTION:
%Latent heat of vaporization in relation to air temperature

         %USES:
%             tfrz = vegetation.physcon.tfrz;         %use clm_varcon                 % Freezing point of water (K)
%             mmh2o = vegetation.physcon.mmh2o;                              %use clm_varcon
%             hfus = vegetation.physcon.hfus; %use clm_varcon
%             hvap = vegetation.physcon.hvap; %use clm_varcon

    %ARGUMENTS:
% t = 200;     % Temperature (K)

    %LOCAL VARIABLES:
%     lambda             % Molar latent heat of vaporization (J/mol)
%---------------------------------------------------------------------

lambda = vegetation.physcon.hvap;                                  % Used in CLM (J/kg)
if (t <= vegetation.physcon.tfrz)
    lambda = lambda + vegetation.physcon.hfus;          % Add latent heat of fusion (J/kg)
end
lambda = lambda * vegetation.physcon.mmh2o;                        % Molar latent heat of vaporization (J/mol)

end
