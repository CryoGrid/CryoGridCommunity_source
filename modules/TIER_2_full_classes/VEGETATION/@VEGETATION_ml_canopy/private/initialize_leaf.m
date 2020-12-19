function [ground] = initialize_leaf(ground)

ground.STATVAR.vegetation.leaf.vpd_min = 100;     % Minimum vapor pressure deficit for Medlyn stomatal conductance (Pa)

% Leaf physiology parameters
% --- Vcmax and other parameters (at 25C)
if (ground.STATVAR.vegetation.pftcon.c3psn == 1)
    ground.STATVAR.vegetation.leaf.vcmax25 = 43;  %https://ldas.gsfc.nasa.gov/gldas/data/GLDAS1_vegparam_tbl_clm2l.pdf
    ground.STATVAR.vegetation.leaf.jmax25 = 1.67 * ground.STATVAR.vegetation.leaf.vcmax25;
    ground.STATVAR.vegetation.leaf.kp25_c4 = 0;
    ground.STATVAR.vegetation.leaf.rd25 = 0.015 * ground.STATVAR.vegetation.leaf.vcmax25;
else
    ground.STATVAR.vegetation.leaf.vcmax25 = 40;
    ground.STATVAR.vegetation.leaf.jmax25 = 0;
    ground.STATVAR.vegetation.leaf.kp25_c4 = 0.02 * ground.STATVAR.vegetation.leaf.vcmax25;
    ground.STATVAR.vegetation.leaf.rd25 = 0.025 * ground.STATVAR.vegetation.leaf.vcmax25;
end

% --- Kc, Ko, Cp at 25C
ground.STATVAR.vegetation.leaf.kc25 = 404.9;
ground.STATVAR.vegetation.leaf.ko25 = 278.4;
ground.STATVAR.vegetation.leaf.cp25 = 42.75;
% --- Activation energy
ground.STATVAR.vegetation.leaf.kcha = 79430;
ground.STATVAR.vegetation.leaf.koha = 36380;
ground.STATVAR.vegetation.leaf.cpha = 37830;
ground.STATVAR.vegetation.leaf.rdha = 46390;
ground.STATVAR.vegetation.leaf.vcmaxha = 65330;
ground.STATVAR.vegetation.leaf.jmaxha  = 43540;
% --- High temperature deactivation
% Deactivation energy (J/mol)
ground.STATVAR.vegetation.leaf.rdhd = 150000;
ground.STATVAR.vegetation.leaf.vcmaxhd = 150000;
ground.STATVAR.vegetation.leaf.jmaxhd  = 150000;
% Entropy term (J/mol/K)
ground.STATVAR.vegetation.leaf.rdse = 490;
ground.STATVAR.vegetation.leaf.vcmaxse = 490;
ground.STATVAR.vegetation.leaf.jmaxse  = 490;
% Scaling factors for high temperature inhibition (25 C = 1.0).
% The factor "c" scales the deactivation to a value of 1.0 at 25C.
% fth25 = @(hd, se) 1 + exp((-hd + se*(physcon.tfrz+25)) / (physcon.rgas*(physcon.tfrz+25)));
ground.STATVAR.vegetation.leaf.vcmaxc = fth25_g (ground, ground.STATVAR.vegetation.leaf.vcmaxhd, ground.STATVAR.vegetation.leaf.vcmaxse);
ground.STATVAR.vegetation.leaf.jmaxc  = fth25_g (ground, ground.STATVAR.vegetation.leaf.jmaxhd, ground.STATVAR.vegetation.leaf.jmaxse);
ground.STATVAR.vegetation.leaf.rdc    = fth25_g (ground, ground.STATVAR.vegetation.leaf.rdhd, ground.STATVAR.vegetation.leaf.rdse);
% --- C3 parameters
% Quantum yield of PS II
ground.STATVAR.vegetation.leaf.phi_psii = 0.85;
% Empirical curvature parameter for electron transport rate
ground.STATVAR.vegetation.leaf.theta_j = 0.90;
% Empirical curvature parameter for C3 co-limitation
ground.STATVAR.vegetation.leaf.colim_c3 = 0.98;
% Empirical curvature parameters for C4 co-limitation
ground.STATVAR.vegetation.leaf.colim_c4a = 0.80;
ground.STATVAR.vegetation.leaf.colim_c4b = 0.95;
% --- C4: Quantum yield (mol CO2 / mol photons)
ground.STATVAR.vegetation.leaf.qe_c4 = 0.05;
% --- Stomatal efficiency for optimization (An/E; umol CO2/ mol H2O)
ground.STATVAR.vegetation.leaf.iota = 750;
% --- ground.STATVAR.vegetation.leaf dimension (m)
% ground.STATVAR.vegetation.leaf.dleaf = 0.05;
% --- ground.STATVAR.vegetation.leaf emissivity
ground.STATVAR.vegetation.leaf.emiss = 0.98;
% --- ground.STATVAR.vegetation.leaf reflectance and transmittance: visible and near-infrared wavebands

ivis = ground.STATVAR.vegetation.params.vis; % Array index for visible waveband
inir = ground.STATVAR.vegetation.params.nir; % Array index for near-infrared waveband
isun = ground.STATVAR.vegetation.params.sun; % Array index for sunlit leaf
isha = ground.STATVAR.vegetation.params.sha; % Array index for shaded leaf

% all values from GLDAS1_vegparameters.pdf
ground.STATVAR.vegetation.leaf.rhol(ivis) = 0.07; %0.10;
ground.STATVAR.vegetation.leaf.taul(ivis) = 0.05; %0.10;
ground.STATVAR.vegetation.leaf.rhol(inir) = 0.35; %0.40;

%%%%%%%%%%%%%%%%%%%%%%%%%%% 20200909 Simone
ground.STATVAR.vegetation.leaf.taul(inir) = 0.25; %0.1;  %0.40;

ground.STATVAR.vegetation.leaf.taus(ivis) = 0.001;
ground.STATVAR.vegetation.leaf.taus(inir) = 0.001;
ground.STATVAR.vegetation.leaf.rhos(ivis) = 0.16;
ground.STATVAR.vegetation.leaf.rhos(inir) = 0.39;

% --- Plant hydraulic parameters  %(Bonan et al. (2014) Geosci. Model Dev., 7, 2193–2222)
% Plant capacitance (mmol H2O/m2 ground.STATVAR.vegetation.leaf area/MPa)
ground.STATVAR.vegetation.leaf.capac = 2500;
% Minimum ground.STATVAR.vegetation.leaf water potential (MPa)
ground.STATVAR.vegetation.leaf.minlwp = -1.5; %-2; %CHANGED SEBASTIAN
% Stem (xylem-to-ground.STATVAR.vegetation.leaf) hydraulic conductance (mmol H2O/m2 ground.STATVAR.vegetation.leaf area/s/Mpa)
ground.STATVAR.vegetation.leaf.gplant = 4.; %0.0004; %4;
end
