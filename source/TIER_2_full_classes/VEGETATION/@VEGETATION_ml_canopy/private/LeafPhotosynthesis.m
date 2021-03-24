function [vegetation] = LeafPhotosynthesis (vegetation, p, ic, il)
% DESCRIPTION:
% Leaf photosynthesis and stomatal conductance

% %LOCAL VARIABLES:
%   qabs                % PAR utilized by PS II (umol photons/m2/s)
%   desat               % Derivative of saturation vapor pressure (Pa/K)
%   gs_err              % gs for error check
%   an_err              % An for error check
%   aquad               % Terms for quadratic equations
%   bquad               % Terms for quadratic equations
%   cquad               % Terms for quadratic equations
%   r1
%   r2                  % Roots of quadratic equation
%   ci0
%   ci1                 % Initial estimates for Ci
%   tol                 % Convergence tolerance for Ci (mol/mol)
%---------------------------------------------------------------------

% *** Input ***
% c3psn  = vegetation.pftcon.c3psn;     % Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant
% g1opt  = vegetation.pftcon.g1opt;     % Ball-Berry slope of conductance-photosynthesis relationship, unstressed
% g0opt  = vegetation.pftcon.g0opt;     % Ball-Berry minimum leaf conductance, unstressed (mol H2O/m2/s)

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%Photosynthesis DONE
%1.2 Leaf photosynthesis and stomatal conductance
%Input
% vegetation.canopy.dpai  = vegetation.canopy.dpai;                     % Layer plant area index (m2/m2)
% btran = vegetation.mlcanopyinst.btran;
% vcmax25 = vegetation.leaf.vcmax25;
% jmax25 = vegetation.leaf.jmax25;
% rd25 = vegetation.leaf.rd25;
% kp25 = vegetation.leaf.kp25;
% eair  = vegetation.mlcanopyinst.eair;                      % Vapor pressure profile (Pa)
% cair  = vegetation.mlcanopyinst.cair;                      % Atmospheric CO2 profile (umol/mol)
% tleaf  = vegetation.mlcanopyinst.tleaf;                    % Leaf temperature (K)
% gbv  = vegetation.mlcanopyinst.gbv;                        % Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
% gbc  = vegetation.mlcanopyinst.gbc;                        % Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
% apar = vegetation.flux.apar;

%Output
%   rd = vegetation.mlcanopyinst.rd;
%   ci = vegetation.mlcanopyinst.ci;
%   hs = vegetation.mlcanopyinst.hs;
%   vpd = vegetation.mlcanopyinst.vpd;

%   %Output from calls to CiFunc or CiFuncGs
%   ac = vegetation.mlcanopyinst.ac;
%   aj = vegetation.mlcanopyinst.aj;
%   ap = vegetation.mlcanopyinst.ap;
%   ag = vegetation.mlcanopyinst.ag;
%   an = vegetation.mlcanopyinst.an;
%   cs = vegetation.mlcanopyinst.cs;
%   gs  = 1;                            % Leaf stomatal conductance (mol H2O/m2 leaf/s)
%   vcmax = vegetation.mlcanopyinst.vcmax;
%   rdleaf = vegetation.mlcanopyinst.rdleaf;
%   jmax = vegetation.mlcanopyinst.jmax;


if (vegetation.canopy.dpai(p,ic) > 0)
    
    %------------------------------------------------------------------
    % Adjust photosynthetic and conductance parameters for temperature
    % and soil water
    %------------------------------------------------------------------
    
    % C3 temperature response
    % ft function
    % how it was: cp     = cp25 .* ft(vegetation, tleaf(p,ic,il), cpha);
    % ERROR: Undefined operator '.*' for input arguments of type 'struct'.
    % So funct result calculated seperately, than multiplied in next
    % step
    %     kcha = vegetation.leaf.kcha;
    %     koha = vegetation.leaf.koha;
    %     cpha = vegetation.leaf.cpha;
    %     vcmaxha = vegetation.leaf.vcmaxha;
    %     vcmaxhd = vegetation.leaf.vcmaxhd;
    %     vcmaxse = vegetation.leaf.vcmaxse;
    %     vcmaxc = vegetation.leaf.vcmaxc;
    %     jmaxha = vegetation.leaf.jmaxha;
    %     jmaxhd = vegetation.leaf.jmaxhd;
    %     jmaxse = vegetation.leaf.jmaxse;
    %     jmaxc = vegetation.leaf.jmaxc;
    %     rdha = vegetation.leaf.rdha;
    %     rdhd = vegetation.leaf.rdhd;
    %     rdse = vegetation.leaf.rdse;
    %     rdc = vegetation.leaf.rdc;
    %     tleaf = vegetation.mlcanopyinst.tleaf(p,ic,il);
    
    [pstr] = ft(vegetation, vegetation.mlcanopyinst.tleaf(p,ic,il), vegetation.leaf.kcha);
    vegetation.leaf.kc     = vegetation.leaf.kc25 .* pstr;
    
    [pstr] = ft(vegetation,vegetation.mlcanopyinst.tleaf(p,ic,il), vegetation.leaf.koha);
    vegetation.leaf.ko     = vegetation.leaf.ko25 .* pstr;
    
    [pstr] = ft(vegetation,vegetation.mlcanopyinst.tleaf(p,ic,il), vegetation.leaf.cpha);
    vegetation.leaf.cp     = vegetation.leaf.cp25 .* pstr;
    
    [pstr] = ft(vegetation,vegetation.mlcanopyinst.tleaf(p,ic,il), vegetation.leaf.vcmaxha);
    [ftv] = fth(vegetation, vegetation.mlcanopyinst.tleaf(p,ic,il), vegetation.leaf.vcmaxhd, vegetation.leaf.vcmaxse, vegetation.leaf.vcmaxc);
    vegetation.leaf.vcmax  = vegetation.leaf.vcmax25(p,ic) .* pstr .* ftv;
    
    [pstr] = ft(vegetation,vegetation.mlcanopyinst.tleaf(p,ic,il), vegetation.leaf.jmaxha);
    [ftv] = fth(vegetation, vegetation.mlcanopyinst.tleaf(p,ic,il), vegetation.leaf.jmaxhd, vegetation.leaf.jmaxse, vegetation.leaf.jmaxc);
    vegetation.leaf.jmax = vegetation.leaf.jmax25(p,ic) .* pstr  .* ftv;
    
    [pstr] = ft(vegetation,vegetation.mlcanopyinst.tleaf(p,ic,il), vegetation.leaf.rdha);
    [ftv] = fth(vegetation,vegetation.mlcanopyinst.tleaf(p,ic,il), vegetation.leaf.rdhd, vegetation.leaf.rdse, vegetation.leaf.rdc);
    vegetation.leaf.rdleaf = vegetation.leaf.rd25(p,ic) .* pstr  .* ftv;
    
    % C4 temperature response
    
    if (round(vegetation.pftcon.c3psn(p)) == 0)  %NINT() -> rounding to the nearest whole number
        vegetation.leaf.vcmax  = vegetation.leaf.vcmax25(p,ic) .* 2.0^((vegetation.mlcanopyinst.tleaf(p,ic,il)-(vegetation.physcon.tfrz+25)) / 10);
        vegetation.leaf.vcmax  = vegetation.leaf.vcmax / ( 1.  + exp(0.2 .*((vegetation.physcon.tfrz+15. )-vegetation.mlcanopyinst.tleaf(p,ic,il))));
        vegetation.leaf.vcmax  = vegetation.leaf.vcmax / ( 1.  + exp(0.3 .*(vegetation.mlcanopyinst.tleaf(p,ic,il)-(vegetation.physcon.tfrz+40))));
        vegetation.leaf.rdleaf = vegetation.leaf.rd25(p,ic) .* 2.0^((vegetation.mlcanopyinst.tleaf(p,ic,il)-(vegetation.physcon.tfrz+25)) / 10);
        vegetation.leaf.rdleaf = vegetation.leaf.rdleaf / ( 1.  + exp(1.3 .*(vegetation.mlcanopyinst.tleaf(p,ic,il)-(vegetation.physcon.tfrz+55))));
    end
    
    vegetation.leaf.kp = vegetation.leaf.kp25(p,ic) .* 2.0.^((vegetation.mlcanopyinst.tleaf(p,ic,il)-(vegetation.physcon.tfrz+25)) / 10);
    % tleaf and leaf.kp25(p,ic)
    
    % Soil water
    
    %     if (vegetation.physcon.gstyp == 0)
    %         vegetation.leaf.vcmax = vegetation.leaf.vcmax .* vegetation.mlcanopyinst.btran(p);
    %         vegetation.mlcanopyinst.g0 = vegetation.pftcon.g0opt(p);
    %     end
   
%         if (gstyp == 1)
%           vegetation.leaf.vcmax = vegetation.leaf.vcmax * vegetation.mlcanopyinst.btran(p);
%           vegetation.mlcanopyinst.g0 = max(vegetation.pftcon.g0opt(p) * vegetation.mlcanopyinst.btran(p), 1.e-06);
%         end 
%        
%         vegetation.leaf.g1 = vegetation.pftcon.g1opt(p);
    
    % btran nicht festsetzen sondern btran skalieren: p.210 or 211, besser p.210
    
    %SEBAS: use constant value of 0.8 set in initialize_mlcanopyinst.m again,
    %probably not used at all, as gstyp is set to 2
    
%     vegetation.mlcanopyinst.btran = 0;
%     
%     if vegetation.soilvar.h2osoi_vol >= 0.35
%         vegetation.mlcanopyinst.btran(p) = 1;
%     elseif vegetation.soilvar.h2osoi_vol <= 0.1 %wilting point
%         vegetation.mlcanopyinst.btran(p) = 0;
%     else
%         vegetation.mlcanopyinst.btran(p) = (vegetation.soilvar.h2osoi_vol(1)-0.1)/(0.35-0.1);
%     end
%         
%     vegetation.mlcanopyinst.btran(p) = max(vegetation.mlcanopyinst.btran(p),1.e-02); %kleiner 1e-02 geht nicht

   %END CHANGE SEBAS
    
    if (vegetation.params.gstyp == 1)
        vegetation.leaf.vcmax = vegetation.leaf.vcmax .* vegetation.mlcanopyinst.btran(p);
        vegetation.mlcanopyinst.g0 = max(vegetation.pftcon.g0opt(p) .* vegetation.mlcanopyinst.btran(p), 1.e-06);
    end
    
    vegetation.leaf.g1 = vegetation.pftcon.g1opt(p);  %g1opt: Ball-Berry slope of conductance-photosynthesis relationship, unstressed
    
    % Save leaf respiration
    vegetation.canopy.rd(p,ic,il) = vegetation.leaf.rdleaf;
    
    %------------------------------------------------------------------
    % Saturation vapor pressure at leaf temperature
    %------------------------------------------------------------------
    
    [esat] = Satvap (vegetation.mlcanopyinst.tleaf(p,ic,il));
    
    vegetation.flux.esat = esat;
    
    %%%%%%%%%% ML & SS Constrain eair >= 0.05.*esat[tleaf] so that solution does not blow up. This
    % ensures that hs does not go to zero. Also eair <= esat[tleaf] so that hs <= 1.
    
    %      ceair = min( max(eair(p,ic), 0.05 .*esat), esat )
    vegetation.leaf.ceair = min(max(vegetation.mlcanopyinst.eair(p,ic), 0.20 .*esat), esat );
    
    %------------------------------------------------------------------
    % Electron transport rate for C3 plants
    %------------------------------------------------------------------
    
    qabs = 0.5  .* vegetation.leaf.phi_psii .* vegetation.flux.apar(p,ic,il);
    aquad = vegetation.leaf.theta_j;
    bquad = -(qabs + vegetation.leaf.jmax);
    cquad = qabs .* vegetation.leaf.jmax;
    [r1,r2] = quadratic(aquad, bquad, cquad);
    vegetation.mlcanopyinst.je = min(r1,r2);
    
    %------------------------------------------------------------------
    % Ci calculation
    %------------------------------------------------------------------
    
    if (vegetation.params.gstyp <= 1)
        
        % Initial estimates for Ci
        
        if (vegetation.pftcon.c3psn(p) == 1)  %NINT() -> rounding to the nearest whole number
            ci0 = 0.7  .* vegetation.mlcanopyinst.cair(p,ic);
        else
            ci0 = 0.4  .* vegetation.mlcanopyinst.cair(p,ic);
        end
        ci1 = ci0 .* 0.99 ;
        
        % Solve for Ci: Use CiFunc to iterate photosynthesis calculations
        % until the change in Ci is < tol. Ci has units umol/mol
        tol = 0.1;
        func_name = 'CiFunc';
        [vegetation, root] = hybrid_root(func_name, vegetation, p, ic, il, ci0, ci1, tol);
        vegetation.mlcanopyinst.ci(p,ic,il) = root;
        
    elseif (vegetation.params.gstyp == 2)
        
        % Calculate photosynthesis for a specified stomatal conductance
     
        [vegetation] = CiFuncGs(vegetation, p,ic,il);
        vegetation.mlcanopyinst.ci(p,ic,il) = vegetation.mlcanopyinst.ci_val;
    end
    
    
    %------------------------------------------------------------------
    % Make sure iterative solution is correct
    %------------------------------------------------------------------
    
    if (vegetation.mlcanopyinst.gs(p,ic,il) < 0. )
        error ('ERROR: LeafPhotosynthesisMod: negative stomatal conductance')
    end
    
    % Compare with Ball-Berry model: vegetation.mlcanopyinst.gs = g1 .* An .* hs/cs + vegetation.mlcanopyinst.g0
    % Use hs calculated with ceair
    
    if (vegetation.params.gstyp == 1)
        vegetation.mlcanopyinst.hs(p,ic,il) = (vegetation.mlcanopyinst.gbv(p,ic,il).*vegetation.leaf.ceair + vegetation.mlcanopyinst.gs(p,ic,il).*esat) ./ ((vegetation.mlcanopyinst.gbv(p,ic,il)+vegetation.mlcanopyinst.gs(p,ic,il)).*esat);
        gs_err = vegetation.leaf.g1.*max(vegetation.mlcanopyinst.an(p,ic,il), 0.).*vegetation.mlcanopyinst.hs(p,ic,il)./vegetation.mlcanopyinst.cs(p,ic,il) + vegetation.mlcanopyinst.g0; %(soil water)
        if (abs(vegetation.mlcanopyinst.gs(p,ic,il)-gs_err)*1.e06  > 1.e-04 )
            disp (' ERROR: LeafPhotosynthesisMod: failed Ball-Berry error check');
        end
    end
    % % %
    % % %     % Compare with Medlyn model: gs = vegetation.mlcanopyinst.g0 + 1.6 .* (1 + g1 / sqrt(Ds)) .* An / cs
    % % %     % Use Ds calculated with ceair and also (esat - ceair) > vpd_min
    % % %
    % % %     if (vegetation.physcon.gstyp == 0)
    % % %         vegetation.leaf.ceair = min(vegetation.leaf.ceair, esat - vpd_min);
    % % %         vegetation.mlcanopyinst.hs(p,ic,il) = (vegetation.mlcanopyinst.gbv(p,ic,il).*vegetation.leaf.ceair + vegetation.mlcanopyinst.gs(p,ic,il).*esat) ./ ((vegetation.mlcanopyinst.gbv(p,ic,il)+vegetation.mlcanopyinst.gs(p,ic,il)).*esat);
    % % %         vegetation.mlcanopyinst.vpd(p,ic,il) = esat - vegetation.mlcanopyinst.hs(p,ic,il) .* esat;
    % % %         gs_err = vegetation.mlcanopyinst.g0 + 1.6  .* (1.  + g1 / sqrt(vegetation.mlcanopyinst.vpd(p,ic,il).*0.001 )) .* max(an(p,ic,il),0. )./cs(p,ic,il);
    % % %         if (abs(vegetation.mlcanopyinst.gs(p,ic,il)-gs_err)*1.e06  > 1. )
    % % %             error ('ERROR: LeafPhotosynthesisMod: failed Medlyn error check');
    % % %         end
    % % %     end
    
    % Compare with diffusion equation: An = (ca - ci) .* gleaf
    %     an_err = (vegetation.mlcanopyinst.cair(p,ic) - ci(p,ic,il)) ./ (1.  / vegetation.mlcanopyinst.gbc(p,ic,il) + 1.6  ./ vegetation.mlcanopyinst.gs(p,ic,il));
    %vegetation.mlcanopyinst.ci(p,ic,il)
    an_err = (vegetation.mlcanopyinst.cair(p,ic) - vegetation.mlcanopyinst.ci(p,ic,il)) ./ (1.  / vegetation.mlcanopyinst.gbc(p,ic,il) + 1.6  ./ vegetation.mlcanopyinst.gs(p,ic,il));
    if (vegetation.mlcanopyinst.an(p,ic,il) > 0) && (abs(vegetation.mlcanopyinst.an(p,ic,il)-an_err) > 0.01)
        disp ('ERROR: LeafPhotosynthesisMod: failed diffusion error check');
    end
    
    %------------------------------------------------------------------
    % Relative humidity and vapor pressure at leaf surface
    %------------------------------------------------------------------
    
    vegetation.mlcanopyinst.hs(p,ic,il) = (vegetation.mlcanopyinst.gbv(p,ic,il).*vegetation.mlcanopyinst.eair(p,ic) + vegetation.mlcanopyinst.gs(p,ic,il).*esat) ./ ((vegetation.mlcanopyinst.gbv(p,ic,il)+vegetation.mlcanopyinst.gs(p,ic,il)).*esat);
    vegetation.mlcanopyinst.vpd(p,ic,il) = max(esat - vegetation.mlcanopyinst.hs(p,ic,il).*esat, 0.1 );
    
else
    %     % non-leaf layer
    %      vegetation.canopy.rd(p,ic,il) = 0;
    %     if (gstyp <= 1)
    %         [ci(p,ic,il), vegetation] = CiFunc(vegetation, p,ic,il);
    %     elseif (gstyp == 2)
    %        [ci(p,ic,il), vegetation] = CiFuncGs(vegetation, p,ic,il);
    %     end
    
    % non-leaf layer
    
    vegetation.mlcanopyinst.rd(p,ic,il) = 0;
    if (vegetation.params.gstyp <= 1)
        [vegetation, ~] = CiFunc(vegetation, p,ic,il,0.);
    elseif (vegetation.params.gstyp == 2)
        [vegetation] = CiFuncGs(vegetation, p,ic,il);
    end
    
    vegetation.mlcanopyinst.hs(p,ic,il) = 0;
    vegetation.mlcanopyinst.vpd(p,ic,il) = 0;
    
end

end

