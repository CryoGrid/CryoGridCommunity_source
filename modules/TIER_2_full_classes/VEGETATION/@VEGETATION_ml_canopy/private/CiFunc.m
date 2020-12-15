function [vegetation, root] = CiFunc(vegetation,p,ic,il,ci_val)

%(p, ic, il, vegetation, ci_val, clm_varctl, MathToolsMod, pftconMod, PatchType, CanopyFluxesMultilayerType)
%
% %DESCRIPTION:
% Calculate leaf photosynthesis and stomatal conductance for a specified Ci
% (ci_val).    calculate a new Ci from the diffusion equation. This
% function equals zero when Ci has converged to the value that satisfies
% the metabolic, stomatal constraint, and diffusion equations.
%
% %USES:
%     use_colim = clm_varctl.use_colim;
%     gstyp = clm_varctl.gstyp;
%quadratic MathToolsMod;
%pftcon  pftconMod;
%patch PatchType;
%mlcanopy_type CanopyFluxesMultilayerType;

%     % %LOCAL VARIABLES:
%     real(r8) :: aquad,bquad,cquad   % Terms for quadratic equations
%     real(r8) :: r1,r2               % Roots of quadratic equation
%     real(r8) :: ai                  % Intermediate co-limited photosynthesis (umol CO2/m2/s)
%     real(r8) :: gleaf               % Leaf CO2 conductance (mol CO2/m2/s)
%     real(r8) :: cinew               % New value for Ci
%     real(r8) :: ci_dif              % Difference in Ci
%     real(r8) :: term                % Term for Medlyn stomatal conductance
%     real(r8) :: vpd_term            % Vapor pressure deficit for Medlyn stomatal conductance (kPa)

% %ARGUMENTS:
% p = vegetation.mlcanopyinst.p;           % Patch index for CLM g/l/c/p hierarchy
% ic = vegetation.mlcanopyinst.ic;                                 % Canopy layer index
% il = vegetation.mlcanopyinst.il;                                 % Sunlit (1) or shaded (2) leaf index
% ci_val = vegetation.mlcanopyinst.cival; % Input value for Ci (umol/mol)

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

gstyp = vegetation.params.gstyp; %user can select between diffrent models for calculating stomatal conductance (gs) 
use_colim = 1;

% rdleaf = vegetation.canopy.rd;
% vegetation.mlcanopyinst.g0  = 1; % Ball-Berry minimum leaf conductance (mol H2O/m2/s); Is this a dummy value?

if (vegetation.canopy.dpai(p,ic) > 0)   % leaf layer
    
    %------------------------------------------------------------------
    % Metabolic (demand-based) photosynthetic rate
    %------------------------------------------------------------------
    
    if (vegetation.pftcon.c3psn(p) == 1)       
        % C3: Rubisco-limited photosynthesis
        vegetation.mlcanopyinst.ac(p,ic,il) = vegetation.leaf.vcmax .* max(ci_val-vegetation.leaf.cp, 0) ./ (ci_val + vegetation.leaf.kc.*(1. + vegetation.mlcanopyinst.o2ref(p)./vegetation.leaf.ko));
        
        % C3: RuBP-limited photosynthesis
        vegetation.mlcanopyinst.aj(p,ic,il) = vegetation.mlcanopyinst.je .* max(ci_val-vegetation.leaf.cp, 0) ./ (4. .*ci_val + 8.*vegetation.leaf.cp);
        
        % C3: Product-limited photosynthesis
        vegetation.mlcanopyinst.ap(p,ic,il) = 0.0;
        
    else        
        % C4: Rubisco-limited photosynthesis
        vegetation.mlcanopyinst.ac(p,ic,il) = vegetation.leaf.vcmax;
        
        % C4: RuBP-limited photosynthesis
        vegetation.mlcanopyinst.aj(p,ic,il) = qe_c4 .* vegetation.flux.apar(p,ic,il);
        
        % C4: PEP carboxylase-limited (CO2-limited)
        vegetation.mlcanopyinst.ap(p,ic,il) = kp .* max(ci_val, 0);        
    end
    
    % Net photosynthesis as the minimum or co-limited rate    
    if (use_colim == 1)
        
        %First co-limit ac and aj       
        if (vegetation.pftcon.c3psn(p) == 1)
            aquad = 0.98;
        else
            aquad = 0.80;
        end
        bquad = -(vegetation.mlcanopyinst.ac(p,ic,il) + vegetation.mlcanopyinst.aj(p,ic,il));
        cquad = vegetation.mlcanopyinst.ac(p,ic,il) .* vegetation.mlcanopyinst.aj(p,ic,il);
        [r1, r2] = quadratic(aquad, bquad, cquad);
        ai = min(r1,r2);
        
        % Now co-limit using ap, but only for C4 plants        
        if (vegetation.pftcon.c3psn(p) == 0)
            vegetation.mlcanopyinst.ag(p,ic,il) = ai;
        else
            aquad = 0.95;
            bquad = -(ai + vegetation.mlcanopyinst.ap(p,ic,il));
            cquad = ai .* vegetation.mlcanopyinst.ap(p,ic,il);
            [r1, r2] = quadratic(aquad, bquad, cquad);
            vegetation.mlcanopyinst.ag(p,ic,il) = min(r1,r2);
        end
               
    elseif (use_colim == 0)           
        if (vegetation.pftcon.c3psn(p) == 1)
            vegetation.mlcanopyinst.ag(p,ic,il) = min(vegetation.mlcanopyinst.ac(p,ic,il),vegetation.mlcanopyinst.aj(p,ic,il));
        else
            vegetation.mlcanopyinst.ag(p,ic,il) = min(vegetation.mlcanopyinst.ac(p,ic,il),vegetation.mlcanopyinst.aj(p,ic,il),vegetation.mlcanopyinst.ap(p,ic,il));
        end
        
    end
    
    % Prevent photosynthesis from ever being negative 
    vegetation.mlcanopyinst.ac(p,ic,il) = max(vegetation.mlcanopyinst.ac(p,ic,il), 0);
    vegetation.mlcanopyinst.aj(p,ic,il) = max(vegetation.mlcanopyinst.aj(p,ic,il), 0);
    vegetation.mlcanopyinst.ap(p,ic,il) = max(vegetation.mlcanopyinst.ap(p,ic,il), 0);
    vegetation.mlcanopyinst.ag(p,ic,il) = max(vegetation.mlcanopyinst.ag(p,ic,il), 0);
    
    % Net photosynthesis    
    vegetation.mlcanopyinst.an(p,ic,il) = vegetation.mlcanopyinst.ag(p,ic,il) - vegetation.leaf.rdleaf;
    
    %------------------------------------------------------------------
    % CO2 at leaf surface
    %------------------------------------------------------------------
    
    vegetation.mlcanopyinst.cs(p,ic,il) = vegetation.mlcanopyinst.cair(p,ic) - vegetation.mlcanopyinst.an(p,ic,il) ./ vegetation.mlcanopyinst.gbc(p,ic,il);
    vegetation.mlcanopyinst.cs(p,ic,il) = max(vegetation.mlcanopyinst.cs(p,ic,il), 1);
    
    %------------------------------------------------------------------
    % Stomatal constraint function
    %------------------------------------------------------------------
    
    % Ball-Berry stomatal conductance
    % Quadratic gs calculation given An. Valid for An >= 0. With An <= 0, gs = vegetation.mlcanopyinst.g0
    
    if (gstyp == 1)
        if (vegetation.mlcanopyinst.an(p,ic,il) > 0)
            aquad = vegetation.mlcanopyinst.cs(p,ic,il);
            bquad = vegetation.mlcanopyinst.cs(p,ic,il).*(vegetation.mlcanopyinst.gbv(p,ic,il) - vegetation.mlcanopyinst.g0) - vegetation.leaf.g1.*vegetation.mlcanopyinst.an(p,ic,il);
            cquad = -vegetation.mlcanopyinst.gbv(p,ic,il) .* (vegetation.mlcanopyinst.cs(p,ic,il).*vegetation.mlcanopyinst.g0 + vegetation.leaf.g1.*vegetation.mlcanopyinst.an(p,ic,il).*vegetation.leaf.ceair./vegetation.flux.esat);
            [r1, r2] = quadratic(aquad, bquad, cquad);
            vegetation.mlcanopyinst.gs(p,ic,il) = max(r1,r2);
        else
            vegetation.mlcanopyinst.gs(p,ic,il) = vegetation.mlcanopyinst.g0;
        end
    end
    
    % Medlyn stomatal conductance
    % Quadratic gs calculation given An. Valid for An >= 0. With An <= 0, gs = vegetation.mlcanopyinst.g0.
    % Note that vapor pressure deficit is limited to be > vpd_min
    
    %     if (gstyp == 0)
    %         if (vegetation.mlcanopyinst.an(p,ic,il) > 0)
    %             vpd_term = max((esat - ceair), vpd_min) .* 0.001;
    %             term = 1.6  .* vegetation.mlcanopyinst.an(p,ic,il) / vegetation.mlcanopyinst.cs(p,ic,il);
    %             aquad = 1.;
    %             bquad = -(2. .* (vegetation.mlcanopyinst.g0 + term) + (g1 .* term)^2. ./ (vegetation.mlcanopyinst.gbv(p,ic,il) .* vpd_term));
    %             cquad = vegetation.mlcanopyinst.g0 .* vegetation.mlcanopyinst.g0 + (2. .* vegetation.mlcanopyinst.g0 + term .* (1. - g1 .* g1 ./ vpd_term)) .* term;
    %             [r1, r2] = quadratic(aquad, bquad, cquad);
    %             vegetation.mlcanopyinst.gs(p,ic,il) = max(r1,r2);
    %         else
    %             vegetation.mlcanopyinst.gs(p,ic,il) = vegetation.mlcanopyinst.g0;
    %         end
    %     end
    
    %------------------------------------------------------------------
    % Diffusion (supply-based) photosynthetic rate - Calculate Ci
    % from the diffusion rate
    %------------------------------------------------------------------
    
    gleaf = 1./(1./vegetation.mlcanopyinst.gbc(p,ic,il) + 1.6 ./ vegetation.mlcanopyinst.gs(p,ic,il));
    vegetation.mlcanopyinst.cinew = vegetation.mlcanopyinst.cair(p,ic) - vegetation.mlcanopyinst.an(p,ic,il)./gleaf;
    
    %------------------------------------------------------------------
    % CiFunc is the difference between the current Ci and the new Ci
    %------------------------------------------------------------------
    
    %     ci_dif(p,ic,il) = cinew - ci_val;
    %     if (an(p,ic,il) < 0)
    %         ci_dif(p,ic,il) = 0;
    
    vegetation.mlcanopyinst.ci_dif = vegetation.mlcanopyinst.cinew - ci_val;
    if (vegetation.mlcanopyinst.an(p,ic,il) < 0)
        vegetation.mlcanopyinst.ci_dif = 0;
    end
    
else % non-leaf layer
    
    vegetation.mlcanopyinst.ac(p,ic,il) = 0;
    vegetation.mlcanopyinst.aj(p,ic,il) = 0;
    vegetation.mlcanopyinst.ap(p,ic,il) = 0;
    vegetation.mlcanopyinst.ag(p,ic,il) = 0;
    vegetation.mlcanopyinst.an(p,ic,il) = 0;
    vegetation.mlcanopyinst.cs(p,ic,il) = 0;
    vegetation.mlcanopyinst.gs(p,ic,il) = 0;
    vegetation.mlcanopyinst.ci_dif = 0;
    
end

root = vegetation.mlcanopyinst.ci_dif;

% Output
% vegetation.mlcanopyinst.ac = ac;
% vegetation.mlcanopyinst.aj = aj;
% vegetation.mlcanopyinst.ap = ap;
% vegetation.mlcanopyinst.ag = ag;
% vegetation.mlcanopyinst.an = an;
% vegetation.mlcanopyinst.cs = cs;
% vegetation.mlcanopyinst.gs = gs;
%ci_dif

end


  
  
  
  
  
  
  
