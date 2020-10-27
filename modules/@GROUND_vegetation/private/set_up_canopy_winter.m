function ground = set_up_canopy_winter(ground)

for p = 1:1 %vegetation.params.npts number of grid points = 1
    
    ground.STATVAR.vegetation.canopy.lai_winter(p) = ground.PARA.lai_winter(p); %4; %5; %8.; %5.0;                                           % Leaf area index of canopy (m2/m2)
    
    % Plant area index of canopy (m2/m2)
    ground.STATVAR.vegetation.canopy.pai_winter(p) = ground.PARA.lai_winter(p) + ground.PARA.sai; %4.051612734794617; %5.051612734794617;
    
    % Atmospheric forcing reference height (m)
    ground.STATVAR.vegetation.mlcanopyinst.zref(p) = ground.PARA.zref; %10.; %15;
    ground.STATVAR.vegetation.mlcanopyinst.sai = ground.PARA.sai;                                        %Stem area index of canopy (m2/m2)                                 ground.STATVAR.vegetation input variables
    
    % Canopy height (m)
    ground.STATVAR.vegetation.mlcanopyinst.ztop(p) = ground.PARA.ztop; %8.; %20
    
    %canopy_type = 'dense';   % High leaf area index
    %canopy_type = 'dense';  % Low leaf area index
    
    %%%%%%%%%%%%
    
    ground.STATVAR.vegetation.mlcanopyinst.nveg(p) = ground.STATVAR.vegetation.mlcanopyinst.nveg_summer(p);
    ground.STATVAR.vegetation.canopy.ntop(p) = ground.STATVAR.vegetation.canopy.ntop_summer(p);
    
    % Set canopy LAI, layer LAI increment, and number of layers
    lai_inc_winter = ground.STATVAR.vegetation.canopy.lai_winter(p)/ground.STATVAR.vegetation.mlcanopyinst.nveg(p); %0.5;
%     ground.STATVAR.vegetation.mlcanopyinst.nveg(p) = round(ground.STATVAR.vegetation.canopy.lai_winter(p) / lai_inc);          % Number of leaf layers in canopy
    lai_inc = lai_inc_winter; 
    
    % Minimum number of layers for Norman radiation
    if (ground.STATVAR.vegetation.mlcanopyinst.nveg(p) < 9)
        ground.STATVAR.vegetation.mlcanopyinst.nveg(p) = 9;
%         lai_inc = ground.STATVAR.vegetation.canopy.lai_winter(p) / ground.STATVAR.vegetation.mlcanopyinst.nveg(p);
    lai_inc_winter = round(ground.STATVAR.vegetation.canopy.lai_winter(p)/ground.STATVAR.vegetation.mlcanopyinst.nveg(p)); %0.5;
    lai_inc = lai_inc_winter; 
    end
       
    % Set array indices for within canopy layers
%     ground.STATVAR.vegetation.soilvar.nsoi(p) = 1;                                                                % First layer is soil ! Also set in initialize_soil
%     ground.STATVAR.vegetation.canopy.nbot(p) = ground.STATVAR.vegetation.soilvar.nsoi(p) + 1;                                    % Bottom leaf layer
%     ground.STATVAR.vegetation.canopy.ntop(p) = ground.STATVAR.vegetation.canopy.nbot(p) + ground.STATVAR.vegetation.mlcanopyinst.nveg - 1;      % Top leaf layer
    
    % Set LAI of each layer 
    for ic = ground.STATVAR.vegetation.canopy.nbot(p):ground.STATVAR.vegetation.canopy.ntop(p)
        ground.STATVAR.vegetation.canopy.dlai_winter(p,ic) = lai_inc;                                                    % dlai = Layer leaf area index
    end
    
    % Cumulative leaf area index (from canopy top) at mid-layer
    for ic = ground.STATVAR.vegetation.canopy.ntop(p):-1:ground.STATVAR.vegetation.canopy.nbot(p) %
        if (ic == ground.STATVAR.vegetation.canopy.ntop(p))
            ground.STATVAR.vegetation.canopy.sumlai_winter(p,ic) = 0.5 * ground.STATVAR.vegetation.canopy.dlai_winter(p,ic);                                %sumlai = Cumulative leaf area index (m2/m2) [for nlevcan layers]
        else
            ground.STATVAR.vegetation.canopy.sumlai_winter(p,ic) = ground.STATVAR.vegetation.canopy.sumlai_winter(p,ic+1) + ground.STATVAR.vegetation.canopy.dlai_winter(p,ic);
        end
    end
    
    % Calculate heights at layer interfaces (zw). These are the heights
    % for the conductance between two scalar concentrations. They are
    % defined for ic = nsoi (ground) to ic = ntop (top of the canopy).
    
    ic = ground.STATVAR.vegetation.canopy.ntop(p);
    ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic) = ground.STATVAR.vegetation.mlcanopyinst.ztop(p); %ztop(p)
    for ic = ground.STATVAR.vegetation.canopy.ntop(p)-1: -1: ground.STATVAR.vegetation.canopy.nbot(p)
        ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic) = ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic+1) - lai_inc;
    end
    
    ic = ground.STATVAR.vegetation.soilvar.nsoi(p);
    if (ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic) > 1e-10 || ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic) < 0)
        error('zw improperly defined at ground level')
    end
    
    % Now calculate the above-canopy layers and their heights
    dz_to_zref = ground.STATVAR.vegetation.mlcanopyinst.zref(p) - ground.STATVAR.vegetation.mlcanopyinst.ztop(p);
    n_to_zref = round(dz_to_zref / lai_inc);
    lai_inc = dz_to_zref / n_to_zref;
    ground.STATVAR.vegetation.mlcanopyinst.ncan(p) = ground.STATVAR.vegetation.canopy.ntop(p) + n_to_zref;
    
    ic = ground.STATVAR.vegetation.mlcanopyinst.ncan(p);
    ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic) = ground.STATVAR.vegetation.mlcanopyinst.zref(p);
    for ic = ground.STATVAR.vegetation.mlcanopyinst.ncan(p)-1: -1: ground.STATVAR.vegetation.canopy.ntop(p)+1
        ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic) = ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic+1) - lai_inc;                  %zw = Canopy heights at layer interfaces (m)
    end
    
    % Determine heights of the scalar concentration and scalar source
    % (zs). These are physically centered between the conductance points
    % (i.e., in the middle of the layer).
    
    ic = ground.STATVAR.vegetation.soilvar.nsoi(p);
    ground.STATVAR.vegetation.mlcanopyinst.zs_winter(p,ic) = 0;
    for ic = ground.STATVAR.vegetation.canopy.nbot(p):ground.STATVAR.vegetation.mlcanopyinst.ncan(p)
        ground.STATVAR.vegetation.mlcanopyinst.zs_winter(p,ic) = 0.5 * (ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic) + ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic-1)); %zs = Canopy height for scalar concentration and source (m)
    end
    
    % Determine plant area index increment for each layer by numerically
    % integrating the plant area density (beta distribution) between
    % the bottom and top heights for that layer
    
    pbeta = 3.5; %11.5; %3.5;     % Parameter for beta distribution
    qbeta = 2; %3.5;     % Parameter for beta distribution
    
    for ic = ground.STATVAR.vegetation.canopy.nbot(p):ground.STATVAR.vegetation.canopy.ntop(p)
        zl = ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic-1); % zw = canopy heights at layer interfaces
        zu = ground.STATVAR.vegetation.mlcanopyinst.zw_winter(p,ic);
        
        ground.STATVAR.vegetation.canopy.dpai_winter(p,ic) = 0;
        
        % Numerical integration between zl and zu using 100 sublayers 
        num_int = 100; %100
        dz_int = (zu - zl) / num_int;
        for ic_int = 1:num_int
            
            if (ic_int == 1)
                z_int = zl + 0.5 * dz_int;
            else
                z_int = z_int + dz_int;
            end
            
            % beta distribution probability density function
            
            zrel = min(z_int/ground.STATVAR.vegetation.mlcanopyinst.ztop(p), 1);
            beta_pdf = (zrel^(pbeta-1) * (1 - zrel)^(qbeta-1)) / beta(pbeta,qbeta);
            
            % Plant area density (m2/m3)
            
            pad = (ground.STATVAR.vegetation.canopy.pai_winter(p) / ground.STATVAR.vegetation.mlcanopyinst.ztop(p)) * beta_pdf;
            
            % Plant area index (m2/m2)
            
            ground.STATVAR.vegetation.canopy.dpai_winter(p,ic) = ground.STATVAR.vegetation.canopy.dpai_winter(p,ic) + pad * dz_int;
            % Plant area index of lowest ground.STATVAR.vegetation layer set, because the whole rest is always left there which causes issues
%             ground.STATVAR.vegetation.canopy.dpai(p,ground.STATVAR.vegetation.canopy.nbot) = 0.2; %2; %32;
            
        end
    end
    
    % Check to make sure sum of numerical integration matches canopy plant area index
    pai_sum = 0;
    for ic = ground.STATVAR.vegetation.canopy.nbot(p):ground.STATVAR.vegetation.canopy.ntop(p)
        pai_sum = pai_sum + ground.STATVAR.vegetation.canopy.dpai_winter(p,ic);
    end
    
       if (abs(pai_sum - ground.STATVAR.vegetation.canopy.pai_winter(p)) > 1e-02) %1e-06
          error('plant area index error')
       end
    
    % Set layers with small plant area index to zero  
    pai_miss = 0;
    for ic = ground.STATVAR.vegetation.canopy.nbot(p):ground.STATVAR.vegetation.canopy.ntop(p)
        if ground.STATVAR.vegetation.canopy.dpai_winter(p,ic) < 0.01 %here the top leaf layer gets PAI = 0. if < 0.01
            pai_miss = pai_miss + ground.STATVAR.vegetation.canopy.dpai_winter(p,ic);
            ground.STATVAR.vegetation.canopy.dpai_winter(p,ic) = 0;
        end
    end
    
    % Distribute the missing plant area across ground.STATVAR.vegetation layers
    % in proportion to the plant area profile
    if (pai_miss > 0)
        pai_old = pai_sum;
        pai_new = pai_old - pai_miss;
        for ic = ground.STATVAR.vegetation.canopy.nbot(p):ground.STATVAR.vegetation.canopy.ntop(p)
            ground.STATVAR.vegetation.canopy.dpai_winter(p,ic) = ground.STATVAR.vegetation.canopy.dpai_winter(p,ic) + pai_miss * (ground.STATVAR.vegetation.canopy.dpai_winter(p,ic) / pai_new);
        end
    end
    
    % Find the lowest ground.STATVAR.vegetation layer   
    for ic = ground.STATVAR.vegetation.canopy.ntop(p): -1: ground.STATVAR.vegetation.canopy.nbot(p)
        if (ground.STATVAR.vegetation.canopy.dpai_winter(p,ic) > 0)
            ic_bot = ic;
        end
    end
    ground.STATVAR.vegetation.canopy.nbot(p) = ic_bot;
    
    % Zero out non-ground.STATVAR.vegetation layers 
    ic = ground.STATVAR.vegetation.soilvar.nsoi(p);
    ground.STATVAR.vegetation.canopy.dpai_winter(p,ic) = 0;
    
    for ic = ground.STATVAR.vegetation.canopy.ntop(p)+1:ground.STATVAR.vegetation.mlcanopyinst.ncan(p)
        ground.STATVAR.vegetation.canopy.dpai_winter(p,ic) = 0;
    end
    
    ground.STATVAR.vegetation.canopy.lai_winter = ground.STATVAR.vegetation.canopy.pai_winter-ground.STATVAR.vegetation.mlcanopyinst.sai; %4.0;
    ground.STATVAR.vegetation.mlcanopyinst.sai = ground.STATVAR.vegetation.canopy.pai_winter - ground.STATVAR.vegetation.canopy.lai_winter;
    ground.STATVAR.vegetation.canopy.dlai_winter = ground.STATVAR.vegetation.canopy.dpai_winter * ground.STATVAR.vegetation.canopy.lai_winter/ground.STATVAR.vegetation.canopy.pai_winter;
    
end
%vegetation.canopy.dlai(:,:) = 0;
%vegetation.canopy.lai = 0;
%---------------------------------------------------------------------
% Cumulative plant area index (lai+sai)
%---------------------------------------------------------------------
ground.STATVAR.vegetation.canopy.num_exposedvegp = 1;
ground.STATVAR.vegetation.canopy.filter_exposedvegp = 1;

for f = 1:ground.STATVAR.vegetation.canopy.num_exposedvegp
    p = ground.STATVAR.vegetation.canopy.filter_exposedvegp(f);
    
    % Layers above the canopy have no ground.STATVAR.vegetation 
    for ic = ground.STATVAR.vegetation.canopy.ntop(p)+1:ground.STATVAR.vegetation.mlcanopyinst.ncan(p)
        ground.STATVAR.vegetation.mlcanopyinst.sumpai_winter(p,ic) = 0.;
    end
    
    % Fill in canopy layers (at the midpoint), starting from the top
    for ic = ground.STATVAR.vegetation.canopy.ntop(p):-1:ground.STATVAR.vegetation.canopy.nbot(p)
        if (ic == ground.STATVAR.vegetation.canopy.ntop(p))
            ground.STATVAR.vegetation.mlcanopyinst.sumpai_winter(p,ic) = 0.5  * ground.STATVAR.vegetation.canopy.dpai_winter(p,ic);
        else
            ground.STATVAR.vegetation.mlcanopyinst.sumpai_winter(p,ic) = ground.STATVAR.vegetation.mlcanopyinst.sumpai_winter(p,ic+1) + 0.5  * (ground.STATVAR.vegetation.canopy.dpai_winter(p,ic+1) + ground.STATVAR.vegetation.canopy.dpai_winter(p,ic));
        end
    end
    
    
ground.STATVAR.vegetation.flux.albsoib_winter = [0.15,0.15]; % Direct beam albedo of ground (soil)
ground.STATVAR.vegetation.flux.albsoid_winter = [0.15,0.15]; % Diffuse albedo of ground (soil)

    for ic = ground.STATVAR.vegetation.canopy.nbot(p): ground.STATVAR.vegetation.canopy.ntop(p)
            ground.STATVAR.vegetation.mlcanopyinst.lwp(p,ic) = -0.1;
            ground.STATVAR.vegetation.mlcanopyinst.h2ocan(p,ic) = 0.;
    end
    
end
