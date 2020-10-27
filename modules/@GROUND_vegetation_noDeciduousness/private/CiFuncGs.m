function [vegetation] = CiFuncGs(vegetation, p,ic,il)

use_colim = 1;
%(p, ic, il, mlcanopy_inst, vcmax, pftcon, clm_varctl, MathToolsMod, pftconMod, PatchType, CanopyFluxesMultilayerType)
% %DESCRIPTION:
% Calculate leaf photosynthesis for a specified stomatal conductance.
% Then calculate Ci from the diffusion equation.
%
% %USES:
% use_colim = 1; %clm_varctl
% quadratic = 1; %MathToolsMod
% pftcon = 1; %pftconMod
% patch = 1; %PatchType
% mlcanopy_type = 1; %CanopyFluxesMultilayerType
%
% %ARGUMENTS:
%     type(mlcanopy_type), intent(inout) :: mlcanopy_inst
%
% %LOCAL VARIABLES:
%     gleaf               % Leaf CO2 conductance (mol CO2/m2/s)
%     a0
%     e0
%     d0            % Terms for quadratic photosynthesis calculation
%     aquad
%     bquad
%     cquad   % Terms for quadratic equations
% r1 = 1;
% r2 = 1;               % Roots of quadratic equation
%     ai = 1;                  % Intermediate co-limited photosynthesis (umol CO2/m2/s)
%     ci_val              % Calculated value for Ci (umol/mol)
%---------------------------------------------------------------------
% *** Input ***

% %ARGUMENTS:
% p = vegetation.mlcanopyinst.p;           % Patch index for CLM g/l/c/p hierarchy
% ic = vegetation.mlcanopyinst.ic;                                 % Canopy layer index
% il = vegetation.mlcanopyinst.il;                                 % Sunlit (1) or shaded (2) leaf index
% nleaf = vegetation.mlcanopyinst.nleaf;
% ncan = vegetation.mlcanopyinst.ncan;
%---------------------------------------------------------------------
% *** Input ***
% c3psn = vegetation.pftcon.c3psn; % Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant
% dpai = vegetation.canopy.dpai; % Layer plant area index (m2/m2)
% o2ref = vegetation.mlcanopyinst.o2ref; % Atmospheric O2 at reference height (mmol/mol)
% cair = vegetation.mlcanopyinst.cair; % Atmospheric CO2 profile (umol/mol)
% gbc = vegetation.mlcanopyinst.gbc; % Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
% gbv  = vegetation.mlcanopyinst.gbv; %Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
% apar = vegetation.flux.apar; % Leaf absorbed PAR (umol photon/m2 leaf/s)

% *** Output ***
% ac = vegetation.mlcanopyinst.ac; %Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
% aj = vegetation.mlcanopyinst.aj; %Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
% ap = vegetation.mlcanopyinst.ap; %Leaf product-limited (C3), CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
% ag = vegetation.mlcanopyinst.ag; %Leaf gross photosynthesis (umol CO2/m2 leaf/s)
% an = vegetation.mlcanopyinst.an; %Leaf net photosynthesis (umol CO2/m2 leaf/s)
% cs = vegetation.mlcanopyinst.cs; %Leaf surface CO2 (umol/mol)
% gs = vegetation.mlcanopyinst.gs; %Leaf stomatal conductance (mol H2O/m2 leaf/s)

% oi    = 1; % Intercellular O2 concentration (mmol/mol)
% j     = 1; % Electron transport rate (umol/m2/s)
% rdleaf = 1;     % Leaf respiration rate (umol CO2/m2 leaf/s)
% g0  = 1;          % Ball-Berry minimum leaf conductance (mol H2O/m2/s)
% g1  = 1;          % Ball-Berry slope of conductance-photosynthesis relationship
% jmax    = 1;     % Maximum electron transport rate (umol/m2/s)

%---------------------------------------------------------------------
% Calculate leaf photosynthesis for a specified stomatal conductance.
% Then calculate Ci from the diffusion equation.
%
% This routine uses a quadratic equation to solve for net photosynthesis (An).
% A general equation for C3 photosynthesis is:
%
%      a (Ci - Cp)
% An = ----------- - Rd
%        e Ci + d
%
% where:
%
% An = Net leaf photosynthesis (umol CO2/m2/s)
% Rd = Leaf respiration (umol CO2/m2/s)
% Ci = Intercellular CO2 concentration (umol/mol)
% Cp = CO2 compensation point (umol/mol)
%
% Rubisco-limited photosynthesis (Ac)
% a  = Vcmax
% e  = 1
% d  = Kc (1 + Oi/Ko)
%
% RuBP-limited photosynthesis (Aj)
% a = J
% e = 4
% d = 8 Cp
%
% where:
%
% Vcmax = Maximum carboxylation rate (umol/m2/s)
% Kc    = Michaelis-Menten constant for CO2 (umol/mol)
% Ko    = Michaelis-Menten constant for O2 (mmol/mol)
% Oi    = Intercellular O2 concentration (mmol/mol)
% J     = Electron transport rate (umol/m2/s)
%
% Ci is calculated from the diffusion equation:
%
%                   1.4   1.6
% An = (Ca - Ci) / (--- + ---)
%                   gb    gs
%
%            1.4   1.6
% Ci = Ca - (--- + ---) An
%            gb    gs
%
% where:
%
% Ca  = Atmospheric CO2 concentration (umol/mol)
% gb  = Leaf boundary layer conductance (mol H2O/m2/s)
% gs  = Leaf stomatal conductance (mol H2O/m2/s)
% 1.4 = Corrects gb for the diffusivity of CO2 compared with H2O
% 1.6 = Corrects gs for the diffusivity of CO2 compared with H2O
%
% The resulting quadratic equation is: a An**2 + b An + c = 0
%
% A similar approach is used for C4 photosynthesis
%---------------------------------------------------------------------

if (vegetation.canopy.dpai(p,ic) > 0)   % leaf layer
    
    % Leaf conductance: gbc has units mol CO2/m2/s, gs has units mol H2O/m2/s,
    % gleaf has units mol CO2/m2/s
    
    gleaf = 1.  / (1. /vegetation.mlcanopyinst.gbc(p,ic,il) + 1.6 /vegetation.mlcanopyinst.gs(p,ic,il));
    
    %------------------------------------------------------------------
    % Gross assimilation rates
    %------------------------------------------------------------------
    
    if (round(vegetation.pftcon.c3psn(p) == 1))
        
        % C3: Rubisco-limited photosynthesis
        
        a0 = vegetation.leaf.vcmax;
        e0 = 1 ;
        d0 = vegetation.leaf.kc .* (1  + vegetation.mlcanopyinst.o2ref(p) / vegetation.leaf.ko);
        
        aquad = e0 / gleaf;
        bquad = -(e0.*vegetation.mlcanopyinst.cair(p,ic) + d0) - (a0 - e0.* vegetation.leaf.rdleaf) / gleaf;
        cquad = a0 .* (vegetation.mlcanopyinst.cair(p,ic) - vegetation.leaf.cp) - vegetation.leaf.rdleaf .* (e0.*vegetation.mlcanopyinst.cair(p,ic) + d0);
        
        
        [r1, r2] = quadratic (aquad, bquad, cquad);
        vegetation.mlcanopyinst.ac(p,ic,il) = min(r1,r2) + vegetation.leaf.rdleaf;
        
        % C3: RuBP-limited photosynthesis
        
        a0 = vegetation.mlcanopyinst.je;
        e0 = 4. ;
        d0 = 8.  .* vegetation.leaf.cp;
        
        aquad = e0 / gleaf;
        bquad = -(e0.*vegetation.mlcanopyinst.cair(p,ic) + d0) - (a0 - e0.*vegetation.leaf.rdleaf) / gleaf;
        cquad = a0 .* (vegetation.mlcanopyinst.cair(p,ic) - vegetation.leaf.cp) - vegetation.leaf.rdleaf .* (e0.*vegetation.mlcanopyinst.cair(p,ic) + d0);
        
        [r1, r2] = quadratic (aquad, bquad, cquad);
        vegetation.mlcanopyinst.aj(p,ic,il) = min(r1,r2) + vegetation.leaf.rdleaf;
        
        % C3: Product-limited photosynthesis
        
        vegetation.mlcanopyinst.ap(p,ic,il) = 0;
        
    else
        
        % C4: Rubisco-limited photosynthesis
        vegetation.mlcanopyinst.ac(p,ic,il) = vegetation.leaf.vcmax;
        
        % C4: RuBP-limited photosynthesis
        vegetation.mlcanopyinst.aj(p,ic,il) = qe_c4 .* vegetation.flux.apar(p,ic,il);
        
        % C4: PEP carboxylase-limited (CO2-limited)
        vegetation.mlcanopyinst.ap(p,ic,il) = vegetation.leaf.kp .* (vegetation.mlcanopyinst.cair(p,ic) .* gleaf + vegetation.leaf.rdleaf) / (gleaf + vegetation.leaf.kp);
        
    end
    
    %------------------------------------------------------------------
    % Net assimilation as the minimum or co-limited rate
    %------------------------------------------------------------------
    
    if (use_colim == 1)
        
        % First co-limit ac and aj
        
        if (round(vegetation.pftcon.c3psn(p)) == 1) %nint()
            aquad = 0.98 ;
        else
            aquad = 0.80 ;
        end
        bquad = -(vegetation.mlcanopyinst.ac(p,ic,il) + vegetation.mlcanopyinst.aj(p,ic,il));
        cquad = vegetation.mlcanopyinst.ac(p,ic,il) .* vegetation.mlcanopyinst.aj(p,ic,il);
        [r1, r2] = quadratic (aquad, bquad, cquad);
        vegetation.mlcanopyinst.ai = min(r1,r2);
        
        % Now co-limit using ap, but only for C4 plants
        
        if (round(vegetation.pftcon.c3psn(p)) == 1) %nint()
            vegetation.mlcanopyinst.ag(p,ic,il) = vegetation.mlcanopyinst.ai;
        else
            aquad = 0.95 ;
            bquad = -(vegetation.mlcanopyinst.ai + vegetation.mlcanopyinst.ap(p,ic,il));
            cquad = vegetation.mlcanopyinst.ai .* vegetation.mlcanopyinst.ap(p,ic,il);
            [r1, r2] = quadratic (aquad, bquad, cquad);
            vegetation.mlcanopyinst.ag(p,ic,il) = min(r1,r2);
        end
        
%     elseif (~ use_colim)
%         
%         if (round(c3psn(p)) == 1)  %nint()
%             ag(p,ic,il) = min(ac(p,ic,il),aj(p,ic,il));
%         else
%             ag(p,ic,il) = min(ac(p,ic,il),aj(p,ic,il),ap(p,ic,il));
%         end
%         
    end
    
    vegetation.mlcanopyinst.an(p,ic,il) = vegetation.mlcanopyinst.ag(p,ic,il) - vegetation.leaf.rdleaf;
    
    %------------------------------------------------------------------
    % Leaf surface CO2
    %------------------------------------------------------------------
    
    vegetation.mlcanopyinst.cs(p,ic,il) = vegetation.mlcanopyinst.cair(p,ic) - vegetation.mlcanopyinst.an(p,ic,il) ./ vegetation.mlcanopyinst.gbc(p,ic,il);
    
    %------------------------------------------------------------------
    % Intercelluar CO2
    %------------------------------------------------------------------
    
    vegetation.mlcanopyinst.ci_val = vegetation.mlcanopyinst.cair(p,ic) - vegetation.mlcanopyinst.an(p,ic,il) ./ gleaf;
    
else
    
    vegetation.mlcanopyinst.ac(p,ic,il) = 0;
    vegetation.mlcanopyinst.aj(p,ic,il) = 0;
    vegetation.mlcanopyinst.ap(p,ic,il) = 0;
    vegetation.mlcanopyinst.ag(p,ic,il) = 0;
    vegetation.mlcanopyinst.an(p,ic,il) = 0;
    vegetation.mlcanopyinst.cs(p,ic,il) = 0;
    vegetation.mlcanopyinst.ci_val = 0;   
end
% Output
%     vegetation.mlcanopyinst.ac = ac;
%     vegetation.mlcanopyinst.aj = aj;
%     vegetation.mlcanopyinst.ap = ap;
%     vegetation.mlcanopyinst.ag = ag;
%     vegetation.mlcanopyinst.an = an;
%     vegetation.mlcanopyinst.cs = cs;
end
