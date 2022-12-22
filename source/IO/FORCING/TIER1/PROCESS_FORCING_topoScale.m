classdef PROCESS_FORCING_topoScale < FORCING_base
    
    
    methods
        
        function forcing = process_topoScale(forcing, tile)
            disp('applying downscaling with TopoScale')
            era = forcing.TEMP.era;
            dist_lat = abs(tile.PARA.latitude - era.lat);
            dist_lon=abs(tile.PARA.longitude-era.lon);
            [dist_lat, ind_lat] = sort(dist_lat);
            [dist_lon, ind_lon] = sort(dist_lon);
            
            dist_lat=dist_lat(1:2);
            dist_lon=dist_lon(1:2);
            
            ind_lat = ind_lat(1:2);
            ind_lon = ind_lon(1:2);
            weights_lat = 1 - dist_lat./sum(dist_lat);
            weights_lon = 1 - dist_lon./sum(dist_lon);
            
            era_alt=era.Z(ind_lon, ind_lat, :, :);
            era_alt=reshape(era_alt, 4, size(era.Z,3), size(era.Z,4));
            
            era_alt_sl  = era.Zs(ind_lon, ind_lat);
            era_alt_sl = reshape(era_alt_sl, 4,1);
            era_alt_sl = repmat(era_alt_sl,1,1,size(era_alt,3));
            
            era_T=era.T(ind_lon, ind_lat, :, :);
            era_T=reshape(era_T, 4, size(era.T,3), size(era.T,4));
            era_u=era.u(ind_lon, ind_lat, :, :);
            era_u=reshape(era_u, 4, size(era.u,3), size(era.u,4));
            era_v=era.v(ind_lon, ind_lat, :, :);
            era_v=reshape(era_v, 4, size(era.v,3), size(era.v,4));
            era_q=era.q(ind_lon, ind_lat, :, :);
            era_q=double(reshape(era_q, 4, size(era.q,3), size(era.q,4))).*era.q_sf;

            
            era_T_sl  = era.T2(ind_lon, ind_lat, :);
            era_T_sl = reshape(era_T_sl, 4,1, size(era_T_sl,3));
            era_wind_sl = sqrt(single(era.u10(ind_lon, ind_lat, :)).^2 + single(era.v10(ind_lon, ind_lat, :)) .^2);
            era_wind_sl = reshape(era_wind_sl, 4,1, size(era_wind_sl,3));
            era_Lin_sl = double(era.LW(ind_lon, ind_lat, :)).*era.rad_sf;
            era_Lin_sl = reshape(era_Lin_sl, 4,1, size(era_Lin_sl,3));
            era_Sin_sl = double(era.SW(ind_lon, ind_lat, :)).*era.rad_sf;
            era_Sin_sl = reshape(era_Sin_sl, 4,1, size(era_Sin_sl,3));
            era_precip_sl = double(era.P(ind_lon, ind_lat, :)).*era.P_sf;
            era_precip_sl = reshape(era_precip_sl, 4,1, size(era_precip_sl,3));
            % era_TOA_sl = era_Sin_sl.*2.5;
            
            %constants
            R=287.05;  % Gas constant for dry air [JK^-1kg^-1]
            g=9.81; % Acceleration of gravity [ms^-1]
            eps0=0.622; % Ratio of molecular weight of water and dry air [-]
            S0=1370; % Solar constat (total TOA solar irradiance) [Wm^-2] used in ECMWF's IFS
            
            K2C=@(Tk) Tk-273.15;
            q2w=@(q) 0.5.*(1-sqrt(1-4.*q)); % Mixing ratio from specific humidity based on 2nd order Taylor series expansion.
            wp2e=@(w,p) 0.5.*p.*(-1+sqrt(1+4.*w./eps0)); % Vapor pressure from mixing ratio based on 2nd order Taylor series expansion.
            % AERK from Alduchov and Eskridge (1996).
            A1=17.625; B1=243.04; C1=610.94;
            Magnus=@(tc) C1.*exp(A1.*tc./(B1+tc)); % A version of the Magnus formula with the AERK parameters.
            % Note, e=Magnus(tdc) and es=Magnus(tc)
            
            era_p_sl = double(era.ps(ind_lon, ind_lat, :)) .* era.ps_sf;
            era_p_sl = reshape(era_p_sl, 4,1, size(era_p_sl,3));
            era_Td_sl = era.Td2(ind_lon, ind_lat, :);
            era_Td_sl = double(reshape(era_Td_sl, 4,1, size(era_Td_sl,3))) .* era.T_sf;
            %vpsl=Magnus(K2C(era_Td_sl));
            vpsl=Magnus(era_Td_sl);
            wsl=eps0.*vpsl./(era_p_sl-vpsl);
            era_q_sl=wsl./(1+wsl);
            
            weights_lat = repmat(weights_lat', 2, 1, 1,size(era_alt,3));
            weights_lat=reshape(weights_lat, 4, 1 , size(era_alt,3));
            weights_lon = repmat(weights_lon, 1, 2, 1, size(era_alt,3));
            weights_lon=reshape(weights_lon, 4, 1 , size(era_alt,3));
            
            layer_below = int16(era_alt.*0);
            layer_above = int16(era_alt.*0);
            
            %do another one to get the lowermost pl above the orography
            
            for i=2:size(era.Z,3)
                layer_below(:,i,:) = double(era_alt(:,i,:) < tile.PARA.altitude & era_alt(:,i-1,:) >= tile.PARA.altitude & era_alt(:,i,:) > era_alt_sl);
                layer_above(:,i-1,:) = double(era_alt(:,i,:) < tile.PARA.altitude & era_alt(:,i-1,:) >= tile.PARA.altitude & era_alt(:,i-1,:) > era_alt_sl);
            end
            
            distance_Z_above = abs(sum(era_alt .* layer_above ,2) - tile.PARA.altitude) .* double(sum(layer_above,2) > 0);
            distance_Z_below = abs(sum(era_alt .* layer_below ,2) - tile.PARA.altitude) .* double(sum(layer_below,2) > 0);
            
            weights_Z_above = 1-distance_Z_above ./ max(1e-10, distance_Z_above + distance_Z_below);
            weights_Z_below = 1-distance_Z_below ./ max(1e-10, distance_Z_above + distance_Z_below);
            weights_Z_above = repmat(weights_Z_above, 1, size(layer_above,2),1) .* double(layer_above);
            weights_Z_below = repmat(weights_Z_below, 1, size(layer_below,2),1) .* double(layer_below);
            
            T_topoScale = sum(double(era_T) .*  (weights_Z_above + weights_Z_below), 2);
            wind_topoScale = sqrt(sum(double(era_u) .*  (weights_Z_above + weights_Z_below), 2).^2 + sum(double(era_v) .*  (weights_Z_above + weights_Z_below), 2).^2) ;
            q_topoScale = sum(double(era_q) .*  double(weights_Z_above + weights_Z_below), 2);
            
            use_sl = sum(weights_Z_above + weights_Z_below,2) <1-1e-9  | tile.PARA.altitude < era_alt_sl;
            merge_w_sl = sum(weights_Z_below,2) == 0;
            factor = min(1, max(0, distance_Z_above./100));
            
            T_topoScale = double(~merge_w_sl) .* T_topoScale + double(merge_w_sl) .* factor .* T_topoScale + double(merge_w_sl) .* (1-factor) .* double(era_T_sl);
            T_topoScale = double(~use_sl) .* T_topoScale + double(use_sl) .* double(era_T_sl);
            wind_topoScale = double(~merge_w_sl) .* wind_topoScale + double(merge_w_sl) .* factor .* wind_topoScale + double(merge_w_sl) .* (1-factor) .* double(era_wind_sl);
            wind_topoScale = double(~use_sl) .* wind_topoScale + double(use_sl) .* double(era_wind_sl);
            q_topoScale = double(~merge_w_sl) .* q_topoScale + double(merge_w_sl) .* factor .* q_topoScale + double(merge_w_sl) .* (1-factor) .* double(era_q_sl);
            q_topoScale = double(~use_sl) .* q_topoScale + double(use_sl) .* double(era_q_sl);
            
            
            
            T_topoScale = T_topoScale .* weights_lat;
            T_topoScale = (T_topoScale(1:2,:,:) +T_topoScale(3:4,:,:));
            wind_topoScale = wind_topoScale .* weights_lat;
            wind_topoScale = (wind_topoScale(1:2,:,:) +wind_topoScale(3:4,:,:));
            q_topoScale = q_topoScale .* double(weights_lat);
            q_topoScale = (q_topoScale(1:2,:,:) + q_topoScale(3:4,:,:));
            era_Lin_sl = era_Lin_sl .* weights_lat;
            era_Lin_sl = (era_Lin_sl(1:2,:,:) + era_Lin_sl(3:4,:,:));
            era_T_sl2 = double(era_T_sl).*era.T_sf .* weights_lat;
            era_T_sl2 = (era_T_sl2(1:2,:,:) + era_T_sl2(3:4,:,:));
            era_alt_sl2 = double(era_alt_sl) .* weights_lat;
            era_alt_sl2 = (era_alt_sl2(1:2,:,:) + era_alt_sl2(3:4,:,:));
            vp_sl2 = double(vpsl) .* weights_lat;
            vp_sl2 = (vp_sl2(1:2,:,:) + vp_sl2(3:4,:,:));
            era_Sin_sl = era_Sin_sl .* weights_lat;
            era_Sin_sl = (era_Sin_sl(1:2,:,:) + era_Sin_sl(3:4,:,:));
            era_precip_sl = era_precip_sl .* weights_lat;
            era_precip_sl = (era_precip_sl(1:2,:,:) + era_precip_sl(3:4,:,:));
            
            weights_lon = (weights_lon(1:2,:,:) + weights_lon(3:4,:,:))./2;
            T_topoScale = squeeze(sum(T_topoScale .* weights_lon,1)) .* era.T_sf;
            wind_topoScale = squeeze(sum(wind_topoScale .* weights_lon,1)) .* era.wind_sf;
            q_topoScale = squeeze(sum(q_topoScale .* double(weights_lon),1));
            era_Lin_sl = squeeze(sum(era_Lin_sl .* weights_lon,1));
            era_alt_sl2 = squeeze(sum(era_alt_sl2 .* weights_lon,1));
            era_T_sl2 = squeeze(sum(era_T_sl2 .* weights_lon,1));
            vp_sl2 = squeeze(sum(vp_sl2 .* weights_lon,1));
            era_Sin_sl = squeeze(sum(era_Sin_sl .* weights_lon,1));
            era_precip_sl = squeeze(sum(era_precip_sl .* weights_lon,1));
            
            
            
            %pressure
            for i=2:size(era.Z,3)
                layer_below(:,i,:) = double(era_alt(:,i,:) < tile.PARA.altitude & era_alt(:,i-1,:) >= tile.PARA.altitude);
                layer_above(:,i-1,:) = double(era_alt(:,i,:) < tile.PARA.altitude & era_alt(:,i-1,:) >= tile.PARA.altitude);
            end
            
            distance_Z_above = abs(sum(era_alt .* layer_above ,2) - tile.PARA.altitude) .* double(sum(layer_above,2) > 0);
            distance_Z_below = abs(sum(era_alt .* layer_below ,2) - tile.PARA.altitude) .* double(sum(layer_below,2) > 0);
            
            weights_Z_above = 1-distance_Z_above ./ max(1e-10, distance_Z_above + distance_Z_below);
            weights_Z_below = 1-distance_Z_below ./ max(1e-10, distance_Z_above + distance_Z_below);
            weights_Z_above = repmat(weights_Z_above, 1, size(layer_above,2),1) .* double(layer_above);
            weights_Z_below = repmat(weights_Z_below, 1, size(layer_below,2),1) .* double(layer_below);
            
            M = 0.0289644;
            R_gc = 8.3144598;
            p_topoScale = sum( repmat(era.p,4,1,size(era.t,2)) .*  (weights_Z_above + weights_Z_below), 2);
            p_topoScale(p_topoScale==0) = era_p_sl(p_topoScale==0) .* exp(-g.*M.*(tile.PARA.altitude-era_alt_sl(p_topoScale==0))./R_gc./(double(era_T_sl(p_topoScale==0)).*era.T_sf + 273.15));
            p_topoScale = p_topoScale .* weights_lat;
            p_topoScale = (p_topoScale(1:2,:,:) + p_topoScale(3:4,:,:));
            p_topoScale = squeeze(sum(p_topoScale .* weights_lon,1));
            
            %needed for Sin
            era_p_sl = era_p_sl .* weights_lat;
            era_p_sl = (era_p_sl(1:2,:,:) + era_p_sl(3:4,:,:));
            era_p_sl = squeeze(sum(era_p_sl .* weights_lon,1));
            
            
            
            %Long-wave
            sbc=5.67e-8;
            era_Lin_sl(era_Lin_sl<=0) = sbc .*(era_T_sl2(era_Lin_sl<=0) + 273.15).^4;
            wf=q2w(q_topoScale); % Convert to mixing ratio at fine grid.
            vpf=wp2e(wf,p_topoScale);
            %disp(qout)
            % Getting negative qout!
            
            % Use the vapor pressure and temperature to calculate clear sky
            % emssivity at grid and subgrid. [also function]
            x1=0.43; x2=5.7;
            cef=real(0.23+x1.*(vpf./(T_topoScale+273.15)).^(1/x2)); % Obs! got negative vapor pressure-> imaginary number in LW calc
            cec=real(0.23+x1.*(vp_sl2./(era_T_sl2+273.15)).^(1/x2));
            
            % Diagnose the all sky emissivity at grid.
            
            aec=era_Lin_sl./(sbc.*(era_T_sl2+273.15).^4);
            
            % Calculate the "cloud" emissivity at grid, assume this is the same at
            % subgrid.
            deltae=aec-cec;
            
            % Use the former cloud emissivity to compute the all sky emissivity at subgrid.
            aef=cef+deltae;
            Lin_topoScale = aef.*sbc.*(T_topoScale+273.15).^4;
            
            
            
            %Sin
            [solar_azm, solar_zen]=solargeom(forcing, era.t,tile.PARA.latitude,tile.PARA.longitude);
            
            solar_azm=[0.5*(solar_azm(1:end-1)+solar_azm(2:end)) solar_azm(end)]';
            solar_zen=[0.5*(solar_zen(1:end-1)+solar_zen(2:end)) solar_zen(end)]';
            
            
            % Compute downwelling TOA SW irradiance (i.e. the incoming shortwave
            % incident on a horizontal plane at the TOA), by accounting for the
            % solar zentih angle.
            mu0=max(cos(solar_zen),0); % Trunacte negative values.
            % May also want to consider treating values mu0<0 for prominent topography
            % when the horizon  angles are less than 0.
            sunset=mu0<cosd(89);%(mu0==0); % Sunset switch.
            % Note, it would be better to use the true average ((1/tau) integral_t^(t+tau) mu0 dt)
            % but this approximation should be ok.
            era_TOA_sl=S0.*mu0;
            
            kt=era_Sin_sl./max(1e-10, era_TOA_sl); % Clearness index.
            kd=0.952-1.041.*exp(-1.*exp(2.3-4.702.*kt)); % Diffuse fraction.
            kd=max(kd,0);
            
            % Diffuse component.
            Sin_diff_topoScale = kd.*era_Sin_sl.*double(~sunset); % Diffuse shortwave radiation.
            % Direct component
            Sin_dir_topoScale = era_Sin_sl - Sin_diff_topoScale;
            
            Sin_dir_topoScale = era_TOA_sl.*(Sin_dir_topoScale./(max(1e-10, era_TOA_sl))).^(p_topoScale./era_p_sl);
            
            % % Scale direct shortwave using Beer's law (see Aalstad 2019, Appendix A)
            % ka=(g.*mu0./(psl)).*log(SWtoa./SWcdir); % Note, we don't get log(0) due to "if sunset" condition.
            % SWfdir=SWtoa.*exp(-ka.*pout./(g*mu0));
            
            
            
            %Precipitation
            threshold_precip = 0.1;
            %  Apply Liston & Elder (MicroMet) elevation-based precip adjustment
            dZ=tile.PARA.altitude - era_alt_sl2; % m
            dZ=dZ./1e3; % km
            dZ=min(dZ,2); % No larger that 2 km=3.3 adjustment
            dZ=max(dZ,-2);% For symmetry
            adjf=0.27; % Mean adjustment factor (following Fiddes and Gruber, 2014)
            adj=(1+adjf.*dZ)./(1-adjf.*dZ);
            
            
            precip_topoScale = era_precip_sl .*adj.*24; %in mm/day, check if timestep must be taken into account when not using 1h input data
            
            
            forcing.DATA.Tair = double(T_topoScale);
            forcing.DATA.q = double(q_topoScale);
            forcing.DATA.wind = double(wind_topoScale);
            forcing.DATA.Sin_dir = double(Sin_dir_topoScale);
            forcing.DATA.Sin_dif = double(Sin_diff_topoScale);
            forcing.DATA.Sin =  forcing.DATA.Sin_dir + forcing.DATA.Sin_dif;
            forcing.DATA.Lin = double(Lin_topoScale);
            forcing.DATA.p = double(p_topoScale);
            forcing.DATA.precip = double(precip_topoScale);
            forcing.DATA.timeForcing = era.t';

        end
        
     
        function forcing = interpolate_sl(forcing, tile)
            disp('interpolating surface level data')
            era = forcing.TEMP.era;
            if length(double(era.lat))>1 && length(double(era.lon))>1
                single_cell=0;
            else
                single_cell=1;
            end
            if single_cell
                ind_lon=1;
                ind_lat=1;
            else
                dist_lat = abs(tile.PARA.latitude - era.lat);
                dist_lon=abs(tile.PARA.longitude-era.lon);
                [dist_lat, ind_lat] = sort(dist_lat);
                [dist_lon, ind_lon] = sort(dist_lon);
                
                dist_lat=dist_lat(1:2);
                ind_lat = ind_lat(1:2);
                weights_lat = 1 - dist_lat./sum(dist_lat);
                dist_lon=dist_lon(1:2);
                ind_lon = ind_lon(1:2);
                weights_lon = 1 - dist_lon./sum(dist_lon);
            end
            
            era_T_sl  = double(era.T2(ind_lon, ind_lat, :)) .* era.T_sf;
            era_wind_sl = sqrt(double(era.u10(ind_lon, ind_lat, :)).^2 + double(era.v10(ind_lon, ind_lat, :)) .^2) .* era.wind_sf;
            era_Lin_sl = double(era.LW(ind_lon, ind_lat, :)).*era.rad_sf;
            era_Sin_sl = double(era.SW(ind_lon, ind_lat, :)).*era.rad_sf;
            era_precip_sl = double(era.P(ind_lon, ind_lat, :)).*era.P_sf;
            era_p_sl = double(era.ps(ind_lon, ind_lat, :)) .* era.ps_sf;
            era_Td_sl = double(era.Td2(ind_lon, ind_lat, :)).* era.T_sf;
            
            if ~single_cell
                era_T_sl = reshape(era_T_sl, 4,1, size(era_T_sl,3));
                era_wind_sl = reshape(era_wind_sl, 4,1, size(era_wind_sl,3));
                era_Lin_sl = reshape(era_Lin_sl, 4,1, size(era_Lin_sl,3));
                era_Sin_sl = reshape(era_Sin_sl, 4,1, size(era_Sin_sl,3));
                era_precip_sl = reshape(era_precip_sl, 4,1, size(era_precip_sl,3));
                era_p_sl = reshape(era_p_sl, 4,1, size(era_p_sl,3));
                era_Td_sl = (reshape(era_Td_sl, 4,1, size(era_Td_sl,3))) ;
            end

            era_q_sl = 0.622 .* (double(era_T_sl>=0) .* satPresWater(forcing, era_Td_sl+273.15) + double(era_T_sl<0) .* satPresIce(forcing, era_Td_sl+273.15)) ./ era_p_sl;
            
            if ~single_cell
                weights_lat = repmat(weights_lat', 2, 1, 1,size(era_T_sl,3));
                weights_lat=reshape(weights_lat, 4, 1 , size(era_T_sl,3));
                weights_lon = repmat(weights_lon, 1, 2, 1, size(era_T_sl,3));
                weights_lon=reshape(weights_lon, 4, 1 , size(era_T_sl,3));
                
                era_wind_sl = era_wind_sl .* weights_lat;
                era_wind_sl = (era_wind_sl(1:2,:,:) +era_wind_sl(3:4,:,:));
                era_q_sl = era_q_sl .* double(weights_lat);
                era_q_sl = (era_q_sl(1:2,:,:) + era_q_sl(3:4,:,:));
                era_Lin_sl = era_Lin_sl .* weights_lat;
                era_Lin_sl = (era_Lin_sl(1:2,:,:) + era_Lin_sl(3:4,:,:));
                era_T_sl = double(era_T_sl) .* weights_lat;
                era_T_sl = (era_T_sl(1:2,:,:) + era_T_sl(3:4,:,:));
                era_Sin_sl = era_Sin_sl .* weights_lat;
                era_Sin_sl = (era_Sin_sl(1:2,:,:) + era_Sin_sl(3:4,:,:));
                era_precip_sl = era_precip_sl .* weights_lat;
                era_precip_sl = (era_precip_sl(1:2,:,:) + era_precip_sl(3:4,:,:));
                era_p_sl = era_p_sl .* weights_lat;
                era_p_sl = (era_p_sl(1:2,:,:) + era_p_sl(3:4,:,:));
                
                weights_lon = (weights_lon(1:2,:,:) + weights_lon(3:4,:,:))./2;
                era_T_sl = squeeze(sum(era_T_sl .* weights_lon,1));
                era_wind_sl = squeeze(sum(era_wind_sl .* weights_lon,1));
                era_q_sl = squeeze(sum(era_q_sl .* double(weights_lon),1));
                era_Lin_sl = squeeze(sum(era_Lin_sl .* weights_lon,1));
                era_Sin_sl = squeeze(sum(era_Sin_sl .* weights_lon,1));
                era_precip_sl = squeeze(sum(era_precip_sl .* weights_lon,1));
                era_p_sl = squeeze(sum(era_p_sl .* weights_lon,1));
            else
                era_T_sl = squeeze(era_T_sl);
                era_wind_sl = squeeze(era_wind_sl);
                era_q_sl = squeeze(era_q_sl);
                era_Lin_sl = squeeze(era_Lin_sl);
                era_Sin_sl = squeeze(era_Sin_sl);
                era_precip_sl = squeeze(era_precip_sl);
                era_p_sl = squeeze(era_p_sl);
            end
            
            forcing.DATA.Tair = double(era_T_sl);
            forcing.DATA.q = double(era_q_sl);
            forcing.DATA.wind = double(era_wind_sl);
            forcing.DATA.Sin =  double(era_Sin_sl);
            forcing.DATA.Lin = double(era_Lin_sl);
            forcing.DATA.p = double(era_p_sl);
            forcing.DATA.precip = double(era_precip_sl) .*24;  %mm/hour to mm/day
            forcing.DATA.timeForcing = era.t';

        end
        
        
    end
    

end

