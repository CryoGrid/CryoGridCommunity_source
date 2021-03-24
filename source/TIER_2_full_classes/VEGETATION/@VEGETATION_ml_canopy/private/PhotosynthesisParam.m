function [vegetation] = PhotosynthesisParam (vegetation, p)

% DESCRIPTION:
% Leaf-level parameters for photosynthesis model

% %USES:
% use_acclim = vegetation.physcon.use_acclim; %clm_varctl

% LOCAL VARIABLES:
%     sco = 1;          % Relative specificity of rubisco

%---------------------------------------------------------------------
%1.1 PhotosynthesisParam (p, mlcanopy_inst)
% Atmospheric O2 at reference height (mmol/mol)
% Average air temperature for acclimation (K)
%---------------------------------------------------------------------
% kc, ko, cp at 25C: Bernacchi et al (2001) Plant, Cell Environment 24:253-259
% Derive sco from cp with o2=0.209 mol/mol and re-calculate cp to allow
% variation in o2
%---------------------------------------------------------------------

vegetation.leaf.kc25 = 404.9;         % umol/mol
vegetation.leaf.ko25 = 278.4;         % mmol/mol
vegetation.leaf.cp25 = 42.75;         % umol/mol

sco = 0.5  .* 0.209  / (vegetation.leaf.cp25 .* 1.e-06 );   % cp25 (umol/mol) -> (mol/mol);
vegetation.leaf.cp25 = 0.5  .* vegetation.mlcanopyinst.o2ref(p) / sco .* 1000;      % O2 is mmol/mol. Multiply by 1000 for umol/mol;

%---------------------------------------------------------------------
% Activation energy:
% Bernacchi et al (2001) Plant, Cell Environment 24:253-259
% Bernacchi et al (2003) Plant, Cell Environment 26:1419-1430

% Acclimation from: Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190
%---------------------------------------------------------------------

vegetation.leaf.kcha    = 79430.;
vegetation.leaf.koha    = 36380.;
vegetation.leaf.cpha    = 37830.;
vegetation.leaf.vcmaxha = 65330.;
vegetation.leaf.jmaxha  = 43540.;
vegetation.leaf.rdha    = 46390.;

% if (use_acclim)
    vegetation.leaf.vcmaxha = 72000.;
    vegetation.leaf.jmaxha  = 50000.;
% end

%---------------------------------------------------------------------
% High temperature deactivation:
% Leuning (2002) Plant, Cell Environment 25:1205-1210
% The factor "c" scales the deactivation to a value of 1.0 at 25C
%
% Acclimation from: Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190
%---------------------------------------------------------------------

vegetation.leaf.vcmaxhd = 150000.;
vegetation.leaf.jmaxhd  = 150000.;
vegetation.leaf.rdhd    = 150000.;

% if (use_acclim)
    vegetation.leaf.vcmaxhd = 200000.;
    vegetation.leaf.jmaxhd  = 200000.;
% end

vegetation.leaf.vcmaxse = 490;
vegetation.leaf.jmaxse  = 490;
vegetation.leaf.rdse    = 490;

% if (use_acclim)
    vegetation.leaf.vcmaxse = 668.39  - 1.07  .* min(max((vegetation.mlcanopyinst.tacclim(p)-vegetation.physcon.tfrz),11.),35.);
    vegetation.leaf.jmaxse  = 659.70  - 0.75  .* min(max((vegetation.mlcanopyinst.tacclim(p)-vegetation.physcon.tfrz),11.),35.);
% end

% fth25 funktion:
[fth25F] = fth25 (vegetation, vegetation.leaf.vcmaxhd, vegetation.leaf.vcmaxse);
vegetation.leaf.vcmaxc = fth25F;

[fth25F] = fth25 (vegetation, vegetation.leaf.jmaxhd, vegetation.leaf.jmaxse);
vegetation.leaf.jmaxc = fth25F;

[fth25F] = fth25 (vegetation, vegetation.leaf.rdhd, vegetation.leaf.rdse);
vegetation.leaf.rdc = fth25F;


%---------------------------------------------------------------------
% Miscellaneous parameters
%---------------------------------------------------------------------

vegetation.leaf.qe_c4 = 0.05 ;
vegetation.leaf.phi_psii = 0.70 ;
%   phi_psii = 0.85
vegetation.leaf.theta_j = 0.90 ;

vegetation.leaf.vpd_min = 100.;
end
