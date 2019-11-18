%-----------------------------------------------------------------------
function [vegetation] = LeafBoundaryLayer (vegetation, p, ic, il)
% (p, ic, il,  pftcon, shr_kind, FORCING)
% %USES:
% tfrz = vegetation.physcon.tfrz;                      % Freezing point of water (K)
% grav = vegetation.physcon.grav;                      % Gravitational acceleration (m/s2)
% rgas = vegetation.physcon.rgas;                      % Universal gas constant (J/K/mol)
% visc0 = vegetation.physcon.visc0;                    % Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
% dh0 = vegetation.physcon.Dh0;                      % Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
% dv0 = vegetation.physcon.Dv0;                      % Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
% dc0 = vegetation.physcon.Dc0;                      % Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)

% % ARGUMENTS:
% nleaf = vegetation.mlcanopyinst.nleaf;
% ncan = vegetation.mlcanopyinst.ncan;
% numrad = vegetation.mlcanopyinst.numrad;
% nlevgrnd = vegetation.mlcanopyinst.nlevgrnd;
% isun = vegetation.params.sun;   %Sunlit leaf
% isha = vegetation.params.sha;   %Shaded leaf
% npts = vegetation.params.npts; %Number of grid points

% p = vegetation.mlcanopyinst.p;           % Patch index for CLM g/l/c/p hierarchy
% ic = 1;%vegetation.mlcanopyinst.ic;                                 % Canopy layer index
% il = 1;%vegetation.mlcanopyinst.il;                                 % Sunlit (1) or shaded (2) leaf index
% dleaf = vegetation.pftcon.dleaf;     % Leaf dimension (m)

%-------------------------------------------------------------------------------
% 1.1 Leaf boundary layer conductance
%                     nlevcan = mlcanopyinst.nlevcan;
%                     nleaf = mlcanopyinst.nleaf;
% vegetation.canopy.dpai(p,ic) = vegetation.canopy.dpai;

% tref(p) = vegetation.mlcanopyinst.tref(p);
% pref(p) = vegetation.mlcanopyinst.pref(p);
% rhomol(p) = vegetation.mlcanopyinst.rhomol(p);
% wind(p,ic) = vegetation.mlcanopyinst.wind(p,ic);

% tleaf(p,ic,il) = vegetation.mlcanopyinst.tleaf(p,ic,il);

% Output
% gbh(p,ic,il) = vegetation.mlcanopyinst.gbh(p,ic,il);
% gbv(p,ic,il) = vegetation.mlcanopyinst.gbv(p,ic,il);
% gbc(p,ic,il) = vegetation.mlcanopyinst.gbc(p,ic,il);

b1 = 1.5; %Empirical correction factor for Nu

if vegetation.canopy.dpai(p,ic) > 0
    
    % Adjust diffusivity for temperature and pressure
    
    fac = 101325.  / vegetation.mlcanopyinst.pref(p) .* (vegetation.mlcanopyinst.tref(p) / vegetation.physcon.tfrz)^1.81 ;
    visc = vegetation.physcon.visc0 .* fac;
    dh = vegetation.physcon.Dh0 .* fac;
    dv = vegetation.physcon.Dv0 .* fac;
    dc = vegetation.physcon.Dc0 .* fac;
    
    % Reynolds number, Prandtl number, Schmidt numbers, and Grashof number
    
    re = vegetation.mlcanopyinst.wind(p,ic) .* vegetation.pftcon.dleaf(p) / visc;
    pr  = visc ./ dh;
    scv = visc ./ dv;
    scc = visc ./ dc;
    gr = vegetation.physcon.grav .* vegetation.pftcon.dleaf(p)^3 .* max(vegetation.mlcanopyinst.tleaf(p,ic,il)-vegetation.mlcanopyinst.tair(p,ic), 0) ./ (vegetation.mlcanopyinst.tair(p,ic) .* visc .* visc);
    
    % Forced convection
    
    % (a) Laminar flow
    
    nu_lam  = b1 .* 0.66   .*  pr^0.33   .* re^0.5  ;
    shv_lam = b1 .* 0.66   .* scv^0.33   .* re^0.5  ;
    shc_lam = b1 .* 0.66   .* scc^0.33   .* re^0.5  ;
    
    % (b) Turbulent flow
    
    nu_turb  = b1 .* 0.036   .*  pr^0.33   .* re^0.8  ;
    shv_turb = b1 .* 0.036   .* scv^0.33   .* re^0.8 ;
    shc_turb = b1 .* 0.036   .* scc^0.33   .* re^0.8  ;
    
    % (c) Choose correct flow regime
    
    nu_forced = max(nu_lam, nu_turb);
    shv_forced = max(shv_lam, shv_turb);
    shc_forced = max(shc_lam, shc_turb);
    
    % Free convection
    
    nu_free  = 0.54   .*  pr^0.25   .* gr^0.25;
    shv_free = 0.54   .* scv^0.25   .* gr^0.25;
    shc_free = 0.54   .* scc^0.25   .* gr^0.25;
    
    % Both forced and free convection regimes occur together
    
    nu = nu_forced + nu_free;
    shv = shv_forced + shv_free;
    shc = shc_forced + shc_free;
    
    % Boundary layer conductances
    
    vegetation.mlcanopyinst.gbh(p,ic,il) = dh .*  nu ./ vegetation.pftcon.dleaf(p);
    vegetation.mlcanopyinst.gbv(p,ic,il) = dv .* shv ./ vegetation.pftcon.dleaf(p);
    vegetation.mlcanopyinst.gbc(p,ic,il) = dc .* shc ./ vegetation.pftcon.dleaf(p);

    % Convert conductance (m/s) to (mol/m2/s)
%     vegetation.mlcanopyinst.rhomol(p) = 40.19999;
    
    vegetation.mlcanopyinst.gbh(p,ic,il) = vegetation.mlcanopyinst.gbh(p,ic,il) .* vegetation.mlcanopyinst.rhomol(p);
    vegetation.mlcanopyinst.gbv(p,ic,il) = vegetation.mlcanopyinst.gbv(p,ic,il) .* vegetation.mlcanopyinst.rhomol(p);
    vegetation.mlcanopyinst.gbc(p,ic,il) = vegetation.mlcanopyinst.gbc(p,ic,il) .* vegetation.mlcanopyinst.rhomol(p);
else  % non-leaf layer
    
    vegetation.mlcanopyinst.gbh(p,ic,il) = 0.;
    vegetation.mlcanopyinst.gbv(p,ic,il) = 0.;
    vegetation.mlcanopyinst.gbc(p,ic,il) = 0.;
    
    % Output
    %         vegetation.mlcanopyinst.gbh = gbh(p,ic,il);
    %         vegetation.mlcanopyinst.gbv = gbv(p,ic,il);
    %         vegetation.mlcanopyinst.gbc = gbc(p,ic,il);
end
end
