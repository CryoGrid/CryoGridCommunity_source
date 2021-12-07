function [T_out, Water_out, WaterIce_out, SnowDepth_out] = CryoGridImplicit(FORCING,GRID,STRAT,PARA,STATVAR,TEMP,SEB)
tic
dt=60*60*24; %-> CONST
t_span = FORCING.t_span; %-> FORCING, full time series, change later in the code
Tair = FORCING.Tair;
snowfall = FORCING.snowfall;
rainfall = FORCING.rainfall;


dxp = GRID.dxp; %-> p: midpoint, delta midpoints , array N-1 x number of tiles
Zp = GRID.Zp; %-> midpoints, cell, array, absolute value of z
Zn = GRID.Zn; % -> upper boundary grid cells, array, absolute value of z
Zs = GRID.Zs; % -> lower boundary all grid cells, array, absolute value of z
dxn = GRID.dxn;  % -> layerThick, N
dxs = GRID.dxs; % -> layerThick, still unclear
kn = GRID.kn;   % -> conductivity through upper boundary
ks = GRID.ks;  % -> conductivity through lower boundary, change in every time step
dxo = GRID.dxo;   % -> lateral distance to the next neighboring tile to the right, used to compute lateral fluxes; array, but stays constant for every cell
An = GRID.An; % -> area of upper interface
As = GRID.As; % -> area of upper interface
Ao = GRID.Ao; % -> area to the right/east, layerThick x contact_length
Vp = GRID.Vp; % -> grid cell volume

cp = TEMP.cp; % -> heat capacity at that point,initialization routuine required that used, see HeatCapacity function below
kp = TEMP.kp; %-> conductivity at that point,initialization routuine required that used, see HeatCapacity function below
lat_flux = TEMP.lat_flux; %-> heat flux at that point, set at zero in the beginning
SnowDepth = TEMP.SnowDepth; %-> one value per cer column, start with zero

WaterIce = STRAT.WaterIce; %-> volumetric fraction 0-1
Mineral = STRAT.Mineral;%-> volumetric fraction 0-1
Organic = STRAT.Organic;%-> volumetric fraction 0-1
Water = STRAT.Water;%-> volumetric fraction 0-1

WaterDepth = PARA.WaterDepth; %-> PARA, in meter, for every column
WaterDensity = PARA.WaterDensity; %-> CONST, column array, 1 x number of tiles
SnowDensity = PARA.SnowDensity; %-> like water depth 1 x number of tiles
Qgeo = PARA.Qgeo; %-> like water depth 1 x number of tiles, W/m2
SnowDepthMax = PARA.SnowCoverMax; %-> like water depth 1 x number of tiles, in meter

T = STATVAR.T; %in degree C
H = STATVAR.H; %energy

N=size(T,1);
J=size(T,2);

T_out = nan(N,length(t_span),J);
Water_out = nan(N,length(t_span),J);
WaterIce_out = nan(N,length(t_span),J);
SnowDepth_out = nan(1,length(t_span),J);

for t = 1:length(t_span)
    fprintf( "Day %d:\n", t );      

    %upper boundary temperature [°C] and snow fall rate [m/day] 
    T0 = Tair(t); %function SEB
        
    % snow and rain precipitation [m/day]
    snow_fall = snowfall(t);
    rain_fall = rainfall(t);
   
    dp = 1.0/dt;
    T_old = T;
    H_old = H;  
  
    for j = J:-1:1 %loop backward through tiles tile 1 is matrix
      
        %Prognostic Snow Scheme
        [Water(:,j), WaterIce(:,j), SnowDepth(1,j), idx] = SnowCover2(Water(:,j),WaterIce(:,j),Zs,SnowDepth(1,j),SnowDepthMax(1,j),snow_fall,SnowDensity,WaterDensity);        
        
        % Bucket hydrology scheme
        %evapotranspiration=0;
        soil_up = sum(Zs<=0.0)+1; %upper most ground cell
        if idx>=soil_up % crude check for snow cover --> needs to be improved
            evapotranspiration = SEB.Qe(j) / (PARA.L_lg * PARA.rho_w) * dt; % in [m] [for one time step], lagging one day behind
            [Water(:,j), WaterIce(:,j)] = BucketScheme(Water(:,j),WaterIce(:,j),Mineral(:,j),Organic(:,j),T(:,j),Zs,dxp,dxs,PARA, rain_fall,evapotranspiration,dt);
        end
        
        %update properties and state variables
        cp(:,j) = HeatCapacity(WaterIce(:,j), Water(:,j), Mineral(:,j), Organic(:,j));
        kp(:,j) = ThermalConductivity(WaterIce(:,j), Water(:,j), Mineral(:,j), Organic(:,j)); 
        [Zn, Zs, dxn, dxs, kn, ks] = MakeGrid(Zp,dxp,kp);       
        H(:,j) = Enthalpie(T(:,j), Water(:,j), cp(:,j));
        [T(:,j), fl] = EnthalpieInv(H(:,j), WaterIce(:,j), cp(:,j));     
        H_old(:,j) = H(:,j);
        T_old(:,j) = T(:,j);
                           
        % Surface energy balance to get upper boundary flux
        temp_offset_max = 1.;
        iter_count = 0;      
        SEB.Tsurf(j) = T(idx,j); % update surface temperature with uppermost grid cell
        %while abs(temp_offset_max)>1e-3 && iter_count<1000
            iter_count=iter_count+1;
            Tsurf_old = SEB.Tsurf(j); 
            Porosity = 1-Mineral-Organic;
            SEB = SurfaceEnergyBalance( T, Water, Porosity, Zp, dxp, kp,  SEB, PARA, FORCING, j, idx, t);
            temp_offset_max = SEB.Tsurf(j) - Tsurf_old;
        %end
        T0=SEB.Tsurf(j);
        
        if iter_count<1000
            fprintf("\t SEB for tile #%d converged after %d iterations: Tsurf=%0.2f, Qg=%0.2f\n", j, iter_count, SEB.Tsurf,SEB.Qg );
        else
            fprintf("Did not converge.\n");
        end
        
        %Prognostic Lake Scheme
        %get the last cell that consists of 100% liquid water
        SolidMatrix = Mineral(:,j) + Organic(:,j);
        lic_lake = Water(idx,j)/(WaterIce(idx,j)+SolidMatrix(idx));
        while lic_lake == 1.0
            %air = 1.0 - WaterIce(idx,j) - SolidMatrix(idx,j);
            lic_lake = Water(idx,j)/(WaterIce(idx,j)+SolidMatrix(idx));
            if lic_lake == 1.0
                idx = idx+1;
            end
        end
            
        %set index of upper boundary
        ubc_idx = idx;

        %Main Iteration Loop
        temp_offset_max = 1;
        iter_count = 1;      
        while (temp_offset_max>=1e-4 || iter_count<2) && iter_count<=1000
            if iter_count>999
                disp('warning: did not reach convergence! Max Temperature offset:')
                disp(temp_offset_max)
            end
            iter_count=iter_count+1;
            
            %implicit update heat capacity and thermal conductivity
            cp(:,j) = HeatCapacity(WaterIce(:,j), Water(:,j), Mineral(:,j), Organic(:,j));
            %kp(:,j) = ThermalConductivity(WaterIce(:,j), Water(:,j), Mineral(:,j), Organic(:,j));
             
            %pre-factores acording to grid cell sizes and thermal conductivities
            kn(ubc_idx) = kp(ubc_idx);%ensures full thermal conductivity at upper boundary
            dxn(ubc_idx) = dxp(ubc_idx)/2.0; %ensures half grid cell at upper boundary
            
            anpn = An(:,j)./Vp(:,j).*kn./dxn;
            anps = As(:,j)./Vp(:,j).*ks./dxs;
            
            %Additional heat fluxes----------------------------------------
            bp = zeros(N,1);
            %flux from upper boundary
            %bp(1) = (An(:,j)./Vp(1,j).*kn(1)./dxn(1)) * T0; %[W/m³]
            %bp(1) = (An(:,j)./Vp(1,j)) * SEB.Qg(j); %[W/m³]
            
            % using surface temperature from SEB:
            bp(1:ubc_idx) = (An(:,j)./Vp(1:ubc_idx,j).*kn(1:ubc_idx)./dxn(1:ubc_idx)) * ( T0 );
            
            % usign ground heat flux from SEB:
            %bp(1:ubc_idx) = (An(:,j)./Vp(1:ubc_idx,j)) * SEB.Qg(j);
            
            anps(1:ubc_idx-1) = 0.0;
            
            
            %flux from lower boundary
            bp(end) = (An(:,j)./Vp(end,j)) * Qgeo; %[W/m³]
            
            
            if j~=1
                ko = (kp(:,1) + kp(:,j))/2.0; %simple mean for testing lat flux
                ap = anpn + anps + Ao(:,j)./Vp(:,j).*ko./dxo(:,j);              
                ap(end) = anpn(end) + Ao(end,j)./Vp(end,j).*ko(end)./dxo(end,j);
                anpn(1:ubc_idx) = 0.0;
                %lateral heat flux (relative [W/m³])
                bp_lat(:,j) = (Ao(:,j)./Vp(:,j).*ko./dxo(:,j)) .* T_old(:,1);
                %update total lateral heat fluxes [W]
                lat_flux(:,j) = (Ao(:,j).*ko./dxo(:,j)) .* (T(:,j)-T_old(:,1));
                %lat_flux(:,j) = (Ao(:,j)./Vp(:,j).*ko./dxo(:,j)) .* (T(:,j)-T_old(:,1));
            else
                %for tile 1 the lateral flux is assumed to be a static
                %external flux reulting from the sum of the lateral heat
                %fluxes from the other tiles.
                ap = anpn + anps;
                ap(end) = anpn(end);   
                
                % modified for Qg
                %ap(1:ubc_idx) = anps(ubc_idx);
                
                
                anpn(1:ubc_idx) = 0.0; 
                bp_lat(:,j) = sum(lat_flux,2)./Vp(:,j); %[W/m³]                
            end
            bp = bp + bp_lat(:,j);
       
            %--------------------------------------------------------------
            [Hinv, ~] = EnthalpieInv(H(:,j), WaterIce(:,j), cp(:,j)); %previous iter
            [~, dHdT] = Enthalpie(Hinv, Water(:,j), cp(:,j)); %previous iter
      
            Sp = -dp*dHdT;
            Sc = dp*(H_old(:,j) - H(:,j)) - Sp.*Hinv;
            
            %TDMA solver --------------------------------------------------
            alpha = -anpn;
            beta = (ap - Sp);
            gamma = -anps;
            delta = Sc + bp;
            T(:,j) = TDMAsolver(alpha,beta,gamma,delta)'; %this iter
            
            %--------------------------------------------------------------
            %update current state of H
            H(:,j) = H(:,j) + dHdT.*(T(:,j) - Hinv);%this iter            
            [Hinv_check, fl] = EnthalpieInv(H(:,j), WaterIce(:,j), cp(:,j));%this iter       
            Water(:,j) = fl.*WaterIce(:,j); %this iter

            R = Hinv_check; %this iter
            B = T(:,j); %this iter
            
            temp_offset = abs(R-B);
            temp_offset_max = max(temp_offset(ubc_idx:end-1,:)); 
            T(:,j) = Hinv_check;
        end
        
        T_out(:,t,j) = T(:,j);
        Water_out(:,t,j) = Water(:,j);
        WaterIce_out(:,t,j) = WaterIce(:,j);
        SnowDepth_out(:,t,j) = SnowDepth(:,j);  
    end
    
    pause(0.01)
    subplot(1,2,1);
    plot(T(:,1),-Zp,T(:,1),-Zp)
    axis([-45 40 -2 0.5])
    grid('on');
    ylabel('depth [m]')
    xlabel('temperature [C]')
    
    pause(0.01)
    subplot(1,2,2);
    plot(Water,-Zp);
    plot(WaterIce,-Zp);
    axis([0 1 -2 0.5])
    grid('on');
    ylabel('depth [m]')
    xlabel('water content [VolFrac]')
end
toc
end



