function [vegetation] = GoudriaanRadiation (omega, kb, kbm, kdm, albcanb, albcand, vegetation)

%     % 
%     % % DESCRIPTION:
%     % Compute solar radiation transfer through canopy using Goudriaan (1977)
%     % as described by Goudriaan and van Laar (1994)
%     % 
%     % % USES:
%     use clm_varpar, only : numrad, isun, isha
%     use clm_varctl, only : iulog
%     % 
%     % % ARGUMENTS:
%     implicit none
%     type(bounds_type), intent(in) :: bounds
%     integer,  intent(in) :: num_exposedvegp                            % Number of non-snow-covered veg points in CLM patch filter
%     integer,  intent(in) :: filter_exposedvegp(:)                      % CLM patch filter for non-snow-covered vegetation
%     real(r8), intent(in) :: omega(bounds%begp:bounds%endp,1:numrad)    % Leaf/stem scattering coefficient
%     real(r8), intent(in) :: kb(bounds%begp:bounds%endp)                % Direct beam extinction coefficient
%     real(r8), intent(in) :: kbm(bounds%begp:bounds%endp,1:numrad)      % kb adjusted for scattering
%     real(r8), intent(in) :: kdm(bounds%begp:bounds%endp,1:numrad)      % kd adjusted for scattering
%     real(r8), intent(in) :: albcanb(bounds%begp:bounds%endp,1:numrad)  % Direct beam albefor above canopy
%     real(r8), intent(in) :: albcand(bounds%begp:bounds%endp,1:numrad)  % Diffuse albefor above canopy
%     type(mlcanopy_type), intent(inout) :: mlcanopy_inst
%     % 
%     % % LOCAL VARIABLES:
%     integer  :: f          % Filter index
%     integer  :: p          % Patch index for CLM g/l/c/p hierarchy
%     integer  :: ic         % Aboveground layer index
%     integer  :: ib         % Waveband index
%     real(r8) :: suminc     % Incident radiation for energy conservation check
%     real(r8) :: sumref     % Reflected radiation for energy conservation check
%     real(r8) :: sumabs     % Absorbed radiation for energy conservation check
% 
%     % Fluxes per unit leaf area (W/m2 leaf)
%     real(r8) :: ild        % Absorbed diffuse flux per unit leaf area (W/m2)
%     real(r8) :: ilb        % Absorbed dir beam (total, with scattering) per unit leaf area (W/m2)
%     real(r8) :: ilbb       % Absorbed direct beam (unscattered) per unit leaf area (W/m2)
%     real(r8) :: ilbs       % Absorbed direct beam (scattered) per unit leaf area (W/m2)
%     real(r8) :: ilsun      % Absorbed solar rad (sunlit leaves) per unit sunlit leaf area (W/m2)
%     real(r8) :: ilsha      % Absorbed solar rad (shaded leaves) per unit shaded leaf area (W/m2)
% 
%     % Fluxes per unit ground area (W/m2 ground area)
%     real(r8) :: icsun      % Absorbed solar radiation, sunlit canopy (W/m2)
%     real(r8) :: icsha      % Absorbed solar radiation, shaded canopy (W/m2)
%     % ---------------------------------------------------------------------
% 
%     associate ( &
%                                                  % *** Input ***
%     clump_fac  => pftcon%clump_fac          , &  % Foliage clumping index (-)
%     vegetation.canopy.lai        => mlcanopy_inst%vegetation.canopy.lai         , &  % Leaf area index of canopy (m2/m2)
%     vegetation.mlcanopyinst.sai        => mlcanopy_inst%vegetation.mlcanopyinst.sai         , &  % Stem area index of canopy (m2/m2)
%     vegetation.mlcanopyinst.ncan       => mlcanopy_inst%vegetation.mlcanopyinst.ncan        , &  % Number of aboveground layers
%     vegetation.canopy.nbot       => mlcanopy_inst%nbot        , &  % Index for bottom leaf layer
%     vegetation.canopy.ntop       => mlcanopy_inst%ntop        , &  % Index for top leaf layer
%     vegetation.canopy.dpai       => mlcanopy_inst%dpai        , &  % Layer plant area index (m2/m2)
%     vegetation.mlcanopyinst.sumpai     => mlcanopy_inst%sumpai      , &  % Cumulative plant area index (m2/m2)
%     vegetation.atmos.swskyb     => mlcanopy_inst%swskyb      , &  % Atmospheric direct beam solar radiation (W/m2)
%     vegetation.atmos.swskyd     => mlcanopy_inst%swskyd      , &  % Atmospheric diffuse solar radiation (W/m2)
%     vegetation.flux.fracsun    => mlcanopy_inst%fracsun     , &  % Sunlit fraction of canopy layer
%     vegetation.flux.fracsha    => mlcanopy_inst%fracsha     , &  % Shaded fraction of canopy layer
%                                                  % *** Output ***
%     vegetation.flux.swleaf     => mlcanopy_inst%swleaf      , &  % Leaf absorbed solar radiation (W/m2 leaf)
%     vegetation.flux.swveg      => mlcanopy_inst%swveg       , &  % Absorbed solar radiation, vegetation (W/m2)
%     vegetation.flux.swvegsun   => mlcanopy_inst%swvegsun    , &  % Absorbed solar radiation, sunlit canopy (W/m2)
%     vegetation.flux.swvegsha   => mlcanopy_inst%swvegsha    , &  % Absorbed solar radiation, shaded canopy (W/m2)
%     albcan     => mlcanopy_inst%albcan      , &  % Albefor above canopy
%     vegetation.flux.swsoi      => mlcanopy_inst%swsoi         &  % Absorbed solar radiation, ground (W/m2)
%     )

isun = vegetation.params.sun; % Array index for sunlit leaf
isha = vegetation.params.sha; % Array index for shaded leaf

clump_fac  = 0.7; %1. ;                     % https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010GB003996  % pft, 1 in matlab exercise Bonan (sp_14_03.m)        ; % Foliage clumping index (-)


    % ---------------------------------------------------------------------
    % Multi-layer radiative transfer
    % ---------------------------------------------------------------------
        
    for ib = 1:vegetation.mlcanopyinst.numrad
       for f = 1:vegetation.canopy.num_exposedvegp
          p = vegetation.canopy.filter_exposedvegp(f);

          % Zero out fluxes for all layers

          for ic = vegetation.mlcanopyinst.nbot(p):vegetation.mlcanopyinst.ncan(p)
             vegetation.flux.swleaf(p,ic,isun,ib) = 0. ;
             vegetation.flux.swleaf(p,ic,isha,ib) = 0. ;
          end

          % Calculate fluxes for leaf layers

          for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
         
             % ild - absorbed diffuse flux per unit leaf area at cumulative vegetation.canopy.lai, 
             % average for all leaves (J / m2 leaf / s)

             ild = (1.  - albcand(p,ib)) * vegetation.atmos.swskyd(p,ib) * kdm(p,ib) * clump_fac(p) ...
                 * exp(-kdm(p,ib) * vegetation.mlcanopyinst.sumpai(p,ic) * clump_fac(p));
            
             % ilb - absorbed direct beam flux (total with scattering) per unit leaf area 
             % at cumulative LAI, average for all leaves (J / m2 leaf / s)

             ilb = (1.  - albcanb(p,ib)) * vegetation.atmos.swskyb(p,ib) * kbm(p,ib) * clump_fac(p) * exp(-kbm(p,ib) * vegetation.mlcanopyinst.sumpai(p,ic) * clump_fac(p));
            
             % ilbb - absorbed direct beam flux (unscattered direct component) per unit leaf area 
             % at cumulative LAI, average for all leaves (J / m2 leaf / s)

             ilbb = (1.  - omega(p,ib)) * vegetation.atmos.swskyb(p,ib) * kb(p) * clump_fac(p) * exp(-kb(p) * vegetation.mlcanopyinst.sumpai(p,ic) * clump_fac(p));
            
             % ilbs - absorbed direct beam flux (scattered direct component) per unit leaf area 
             % at cumulative LAI, average for all leaves (J / m2 leaf / s)

             ilbs = ilb - ilbb;
            
             % ilsha - total absorbed flux (shaded leaves) per unit shaded leaf area 
             % at cumulative LAI (J / m2 leaf / s)

             ilsha = ild + ilbs;
            
             % ilsun - total absorbed flux (sunlit leaves) per unit sunlit leaf area 
             % at cumulative LAI (J / m2 leaf / s)

             ilsun = ilsha + kb(p) * (1.  - omega(p,ib)) * vegetation.atmos.swskyb(p,ib);

             % Solar radiation absorbed by sunlit and shaded leaf

             vegetation.flux.swleaf(p,ic,isun,ib) = ilsun;
             vegetation.flux.swleaf(p,ic,isha,ib) = ilsha;

          end
       end
    end

    % ---------------------------------------------------------------------
    % Canopy summation and soil absorption
    % ---------------------------------------------------------------------

    for ib = 1:vegetation.mlcanopyinst.numrad
       for f = 1:vegetation.canopy.num_exposedvegp
          p = vegetation.canopy.filter_exposedvegp(f);

          % Canopy integration, sunlit and shaded leaves

          icsun = 0. ;
          icsha = 0. ;

          for ic = vegetation.canopy.nbot(p), vegetation.canopy.ntop(p)
             icsun = icsun + vegetation.flux.swleaf(p,ic,isun,ib) * vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic);
             icsha = icsha + vegetation.flux.swleaf(p,ic,isha,ib) * vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic);
          end

          % Solar radiation absorbed by vegetation

          vegetation.flux.swveg(p,ib) = icsun + icsha;
          vegetation.flux.swvegsun(p,ib) = icsun;
          vegetation.flux.swvegsha(p,ib) = icsha;

          % Solar radiation absorbed by ground (soil)

          vegetation.flux.swsoi(p,ib) = vegetation.atmos.swskyb(p,ib) * (1.  - albcanb(p,ib)) * exp(-kbm(p,ib)*(vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p))*clump_fac(p)) ...
                      + vegetation.atmos.swskyd(p,ib) * (1.  - albcand(p,ib)) * exp(-kdm(p,ib)*(vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p))*clump_fac(p));

          % Solar radiation reflected by canopy

          suminc = vegetation.atmos.swskyb(p,ib) + vegetation.atmos.swskyd(p,ib);
          sumabs = vegetation.flux.swveg(p,ib) + vegetation.flux.swsoi(p,ib);
          sumref = suminc - sumabs;
%         sumref = vegetation.atmos.swskyb(p,ib) * albcanb(p,ib) + vegetation.atmos.swskyd(p,ib) * albcand(p,ib)

          if (suminc > 0. ) 
             vegetation.flux.albcan(p,ib) = sumref / suminc;
          else
             vegetation.flux.albcan(p,ib) = 0. ;
          end

       end
    end

    end