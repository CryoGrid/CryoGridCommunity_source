
%========================================================================

classdef FORCING_seb_topoScale_CCI < matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        TEMP
    end
    
    
    methods
        
        
        function forcing = provide_PARA(forcing)
            forcing.PARA.ERA_path = [];
            forcing.PARA.filename_SURF_airT = [];
            forcing.PARA.filename_SURF_geopotential = [];
            forcing.PARA.filename_SURF_dewpointT = [];
            forcing.PARA.filename_SURF_solar_radiation = [];
            forcing.PARA.filename_SURF_thermal_radiation = [];
            forcing.PARA.filename_SURF_TOA_radiation = [];
            forcing.PARA.filename_SURF_precip = [];
            forcing.PARA.filename_PLEV_airT = [];
            forcing.PARA.filename_PLEV_geopotential = [];
            forcing.PARA.filename_PLEV_RH = [];
            
            forcing.PARA.start_time = [];
            forcing.PARA.end_time = [];
            
        end
        
        
        
        function forcing = provide_CONST(forcing)
           forcing.CONST.g =[]; 
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end
        
        function forcing = finalize_init(forcing, tile)
            
            forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1), forcing.PARA.end_time(3,1));
            
            test_year = year(forcing.PARA.start_time); %choose a year that exists to read the time-constant data sets
            filename = forcing.PARA.filename_SURF_airT;
            filename(21:24) = num2str(test_year);

            ERA_lat = ncread([forcing.PARA.ERA_path filename], 'latitude');
            ERA_lon = ncread([forcing.PARA.ERA_path filename], 'longitude');
            ERA_lon =[ERA_lon; single(360)]; %add the 0-degree-longitude to the end as 360-degree-longitude
            
            max_target_lat = min(round(max(tile.PARA.latitude) *4) /4 + 0.25, 90) ;  %ERA5 in 0.25 degree resolution
            min_target_lat= max(round(min(tile.PARA.latitude)*4) /4 - 0.25, -90);
%             max_target_lat = min(floor(max(tile.PARA.latitude+0.25) *4) /4 + 0.25, 90) ;  %ERA5 in 0.25 degree resolution
%             min_target_lat= max(floor(min(tile.PARA.latitude)*4) /4, -90);
            
            target_lon = tile.PARA.longitude;
            target_lon(target_lon < 0) = target_lon(target_lon < 0) + 360; %change from -180 to +180 to 0 to 360
            max_target_lon = min(round(max(target_lon)*4) /4 + 0.25, 360);
            min_target_lon = max(round(min(target_lon)*4) / 4 - 0.25, 0);
            
            forcing.TEMP.lat_index_start = find(ERA_lat(:,1) == max_target_lat);  
            forcing.TEMP.lat_index_end = find(ERA_lat(:,1) == min_target_lat);
            forcing.TEMP.lon_index_start = find(ERA_lon(:,1) == min_target_lon);
            forcing.TEMP.lon_index_end = find(ERA_lon(:,1) == max_target_lon);
            
            forcing.TEMP.max_target_lon = max_target_lon;
            forcing.TEMP.target_lon = target_lon;
            
            filename = forcing.PARA.filename_SURF_geopotential;
            filename(19:22) = num2str(test_year);
            
            forcing.TEMP.z_surface = ncread_with_360_degree(forcing, [forcing.PARA.ERA_path filename], 'z', 1, 1) ./ forcing.CONST.g;
            
            [forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh] = meshgrid(ERA_lon(forcing.TEMP.lon_index_start:forcing.TEMP.lon_index_end), ERA_lat(forcing.TEMP.lat_index_start:forcing.TEMP.lat_index_end));
            
            forcing.TEMP.interpolated_ERA_orography = single(interp2(forcing.TEMP.ERA_lon_mesh , forcing.TEMP.ERA_lat_mesh, forcing.TEMP.z_surface', target_lon, tile.PARA.latitude, 'cubic'));

            %load ERA data
            
            forcing.TEMP.ERA_time=[];

            number_of_values = (forcing.PARA.end_time-forcing.PARA.start_time+1).*4;
            forcing.TEMP.ERA_T_downscaled=single(zeros(size(tile.PARA.latitude,1),number_of_values));
            forcing.TEMP.ERA_Lin_downscaled = single(zeros(size(tile.PARA.latitude,1),number_of_values));
            forcing.TEMP.ERA_Sin_downscaled = single(zeros(size(tile.PARA.latitude,1),number_of_values));
            forcing.TEMP.ERA_precip_downcaled=single(zeros(size(tile.PARA.latitude,1),number_of_values));
            
            count = 0;
            for time = forcing.PARA.start_time:forcing.PARA.end_time
                datestr(time)
                current_year = str2num(datestr(time,'yyyy'));
                doy = time - datenum(current_year,1,1) + 1;
                time_index_start = (doy-1) .* 4 + 1;
                
                filename_SURF_precip = forcing.PARA.filename_SURF_precip;
                filename_SURF_precip(26:29) = num2str(current_year);
                
                filename_SURF_airT = forcing.PARA.filename_SURF_airT;
                filename_SURF_airT(21:24) = num2str(current_year);
                
                filename_SURF_dewpointT = forcing.PARA.filename_SURF_dewpointT;
                filename_SURF_dewpointT(30:33) = num2str(current_year);
                
                filename_SURF_solar_radiation = forcing.PARA.filename_SURF_solar_radiation;
                filename_SURF_solar_radiation(40:43) = num2str(current_year);
                
                filename_SURF_thermal_radiation = forcing.PARA.filename_SURF_thermal_radiation;
                filename_SURF_thermal_radiation(42:45) = num2str(current_year);
                
                filename_SURF_TOA_radiation = forcing.PARA.filename_SURF_TOA_radiation;
                filename_SURF_TOA_radiation(35:38) = num2str(current_year);
                
                filename_PLEV_airT = forcing.PARA.filename_PLEV_airT;
                filename_PLEV_airT(18:21) = num2str(current_year);
                
                filename_PLEV_geopotential = forcing.PARA.filename_PLEV_geopotential;
                filename_PLEV_geopotential(19:22) = num2str(current_year);               
                
                filename_PLEV_RH = forcing.PARA.filename_PLEV_RH;
                filename_PLEV_RH(24:27) = num2str(current_year); 

                if exist([forcing.PARA.ERA_path filename_SURF_airT])==2
                    T2m = ncread_with_360_degree(forcing, [forcing.PARA.ERA_path filename_SURF_airT], 't2m', time_index_start, 4);
                    
                    tot_prec = ncread_with_360_degree(forcing, [forcing.PARA.ERA_path filename_SURF_precip], 'tp', time_index_start, 4); 
                    
                    T_700 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_airT], 't', time_index_start, 4, 2);
                                        
                    z_700 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_geopotential], 'z', time_index_start, 4,2) ./ forcing.CONST.g;
                    
                    T_500 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_airT], 't', time_index_start, 4, 3);
                    
                    z_500 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_geopotential], 'z', time_index_start, 4, 3) ./ forcing.CONST.g;
                    
                    T_200 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_airT], 't', time_index_start, 4, 4);
                    
                    z_200 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_geopotential], 'z', time_index_start, 4, 4) ./ forcing.CONST.g;
                    
                    diff700 = z_700 - forcing.TEMP.z_surface;
                    diff500 = z_500 - forcing.TEMP.z_surface;
                    diff200 = z_200 - forcing.TEMP.z_surface;
                    
                    factor700 = double(diff700>500) .* normpdf(diff700, 3000, 800); %this gives a half-way smooth transition, the 200mbar level is not optimal, more than 10km, should be a bit lower, like 300 or 350mbar
                    factor500 = double(diff500>500) .* normpdf(diff500, 3000, 800);
                    factor200 = double(diff200>500) .* normpdf(diff200, 5000, 1600);
                    
                    sum_factors = factor700+factor500+factor200; %normalize
                    factor700 = factor700./sum_factors;
                    factor500 = factor500./sum_factors;
                    factor200 = factor200./sum_factors;
                    
                    lapse_rate = factor700 .* (T_700-T2m) ./ (z_700 - forcing.TEMP.z_surface) + factor500 .* (T_500-T2m) ./ (z_500 - forcing.TEMP.z_surface) + factor200 .* (T_200-T2m) ./ (z_200 - forcing.TEMP.z_surface);

                    T_dew_ERA = ncread_with_360_degree(forcing, [forcing.PARA.ERA_path filename_SURF_dewpointT], 'd2m', time_index_start, 4);
                    
                    Sin_ERA = ncread_with_360_degree(forcing, [forcing.PARA.ERA_path filename_SURF_solar_radiation], 'ssrd', time_index_start, 4) ./ 3600; 
                    
                    Lin_ERA = ncread_with_360_degree(forcing, [forcing.PARA.ERA_path filename_SURF_thermal_radiation], 'strd', time_index_start, 4) ./ 3600; 
                    
                    TOA_ERA = ncread_with_360_degree(forcing, [forcing.PARA.ERA_path filename_SURF_TOA_radiation], 'tisr', time_index_start, 4) ./ 3600;
                                        
                    RH_700 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_RH], 'r', time_index_start, 4, 2);
                    
                    RH_500 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_RH], 'r', time_index_start, 4, 3);
                    
                    RH_200 = ncread_with_360_degree_pl(forcing, [forcing.PARA.ERA_path filename_PLEV_RH], 'r', time_index_start, 4, 4);
                    
                    A1=7.625; B1=243.04; C1=610.94; Tm=273.15;
                    dew2RH=@(tdc,tc) 1e2.*exp((A1.*tdc-(tdc+B1).*(A1.*tc./(B1+tc)))./(tdc+B1)); % Converts from dewpoint temp (C) to RH (%)
                    vpsat=@(tc) C1.*exp((A1.*tc)./(B1+tc)); % Converts from temp (C) to saturation vapor pressure (Pa)
                    RH2vp=@(rh,tc) rh.*vpsat(tc)./1e2;  % Converts from RH (%) to vapor pressure (Pa).h
%                     RHi=dew2RH(Tdi-Tm,Ti-Tm);
%                     vpi=RH2vp(RHi,Ti-Tm);
%                     vp=RH2vp(RH,T-Tm);
                    
                    RH_ERA_surface = dew2RH(T_dew_ERA - Tm, T2m - Tm);
                    
                    lapse_rate_RH = factor700 .* (RH_700-RH_ERA_surface) ./ (z_700 - forcing.TEMP.z_surface) + factor500 .* (RH_500-RH_ERA_surface) ./ (z_500 - forcing.TEMP.z_surface) + factor200 .* (RH_200-RH_ERA_surface) ./ (z_200 - forcing.TEMP.z_surface);
                    
                    
                    
                    for i=1:4
                        %interpolate airT
                        T2m_slice=squeeze(T2m(:,:,i));
                        interpolated_T2m = interp2(forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh, T2m_slice', forcing.TEMP.target_lon, tile.PARA.latitude, 'cubic');
                        lr_slice=squeeze(lapse_rate(:,:,i)); 
                        interpolated_lapse_rate = interp2(forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh, lr_slice', forcing.TEMP.target_lon, tile.PARA.latitude, 'cubic');
                        
                        T_downscaled =  interpolated_T2m + interpolated_lapse_rate .* (tile.PARA.altitude - forcing.TEMP.interpolated_ERA_orography);
                        
                        %interpolate RH
                        RH_ERA_surface_slice = squeeze(RH_ERA_surface(:,:,i));
                        interpolated_RH_surface = interp2(forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh, RH_ERA_surface_slice', forcing.TEMP.target_lon, tile.PARA.latitude, 'cubic');
                        lr_slice=squeeze(lapse_rate_RH(:,:,i));
                        interpolated_lapse_rate = interp2(forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh, lr_slice', forcing.TEMP.target_lon, tile.PARA.latitude, 'cubic');
                        
                        RH_downscaled =  interpolated_RH_surface + interpolated_lapse_rate .* (tile.PARA.altitude - forcing.TEMP.interpolated_ERA_orography);
                        RH_downscaled = min(100,max(0, RH_downscaled));
                        
                        vpi = RH2vp(interpolated_RH_surface, interpolated_T2m - Tm); %vapor pressure
                        vp = RH2vp(RH_downscaled, T_downscaled - Tm);
                        
                        % Longwave routine
                        % Same as TopoSCALE, but without terrain corrections.
                        
                        % Clear sky emissivity (Konzelmann 1994 parametrization)
                        x1=0.43; x2=5.7;
                        epsclear=@(vpress,tk) 0.23+x1.*(vpress./tk).^(1/x2);
                        emcsi=epsclear(vpi,interpolated_T2m);
                        emcs=epsclear(vp,T_downscaled);
                        
                        
                        Lin_ERA_slice = squeeze(Lin_ERA(:,:,i));
                        interpolated_Lin = interp2(forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh, Lin_ERA_slice', forcing.TEMP.target_lon, tile.PARA.latitude, 'cubic');
                        % All-sky emissivity at zi
                        sbc=5.67e-8;
                        emasi = interpolated_Lin./(sbc.*interpolated_T2m.^(4));
                        
                        % Cloud emissivity (assumed same at zi and z)
                        emcl = emasi - emcsi;
                        
                        % All-sky emissivty at z
                        emas=emcs+emcl;
                        
                        % Downwelling longwave at z
                        Lin_downscaled =emas.*sbc.*T_downscaled.^4;

                        
                        % Shortwave routine
                        % Same as TopoSCALE, but without terrain corrections.
                        
                        % Direct-diffuse partitioning.
                        Sin_ERA_slice = squeeze(Sin_ERA(:,:,i));
                        interpolated_Sin = interp2(forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh, Sin_ERA_slice', forcing.TEMP.target_lon, tile.PARA.latitude, 'cubic');
                        TOA_ERA_slice = squeeze(TOA_ERA(:,:,i));
                        interpolated_TOA = interp2(forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh, TOA_ERA_slice', forcing.TEMP.target_lon, tile.PARA.latitude, 'cubic');

                        kt = interpolated_Sin./interpolated_TOA;
                        kt(isnan(kt))=1;
                        a=0.952; b=1.041; c=2.3; d=4.702;
                        % Double exponential relationship for diffuse fraction
                        % from Ruiz-Arias et al. (2010)
                        fdiff=@(ci) a-b.*exp(-1.*exp(c-d.*ci));
                        fd=fdiff(kt);
                        fd=fd.*(fd>0);
                        Sidiff=fd.*interpolated_Sin;
                        Sidir = interpolated_Sin - Sidiff;
                        
                        M=0.02896; % kg mol?1
                        g=9.81;
                        R=8.314;
                        
                        p1_over_p2_interp = (1 + interpolated_lapse_rate .* (tile.PARA.altitude - forcing.TEMP.interpolated_ERA_orography) ./ interpolated_T2m) .^(M.*g./R./-interpolated_lapse_rate);
                         
                        % Elevation scaling based on Beer's law (see Appendix A in my (Kris') thesis)
                        beer=@(p1_over_p2,Stoa,Sdir2) Stoa.*(Sdir2./Stoa).^(p1_over_p2);
                        Sdir=beer(p1_over_p2_interp, interpolated_TOA, Sidir);
                        Sdir(isnan(Sdir))=0; %polar night
                        % Assuming diffuse radiation is independent with height, this gives us
                        % the following incoming shortwave radiation at z:
                        Sin_downscaled  = real(Sdir + Sidiff);
                        Sin_downscaled(Sin_downscaled<0)=0;
                                                
                        
                        forcing.TEMP.ERA_time=[forcing.TEMP.ERA_time time + (i-1).*0.25];
                        forcing.TEMP.ERA_T_downscaled(:,count+i) = T_downscaled;
                        forcing.TEMP.ERA_Lin_downscaled(:,count+i) = Lin_downscaled;
                        forcing.TEMP.ERA_Sin_downscaled(:,count+i) = Sin_downscaled;
                    end

                    %precip
                    for jjj=1:4  %jjj=1:6:19
                        
                        tot_prec_slice = tot_prec(:,:,jjj) .* 24; % ./ forcing.PARA.number_of_values_per_day; % 6;  %accumulate over 24h slice, so unit is m/day
                        %   tot_prec_slice = sum(tot_prec(:,:,jjj:jjj+5),3);  %accumulate over 6h slice
                        interpolated_tot_prec = interp2(forcing.TEMP.ERA_lon_mesh, forcing.TEMP.ERA_lat_mesh, tot_prec_slice', forcing.TEMP.target_lon, tile.PARA.latitude, 'cubic');
                        
                        tot_prec_downscaled = interpolated_tot_prec .* 1.04.^((tile.PARA.altitude - forcing.TEMP.interpolated_ERA_orography)./100);  %5 percent increase per 100m elevation increase -> change to Jaros formula later
                        tot_prec_downscaled(tot_prec_downscaled<0) = 0;
                        %append values to file
                        %forcing.TEMP.ERA_precip_downcaled =[forcing.TEMP.ERA_precip_downcaled tot_prec_downscaled.*1000]; %mm/day
                        
                        forcing.TEMP.ERA_precip_downcaled(:,count+jjj) = tot_prec_downscaled.*1000; %mm/day
                    end
                    
                    
                else
                    forcing.TEMP.ERA_time=[forcing.TEMP.ERA_time time + [0 0.25 0.5 0.75]];
                    
                    forcing.TEMP.ERA_T_downscaled(:,count+1:count+4)= repmat(273.15-10,size(tile.PARA.latitude,1), 4);
                    forcing.TEMP.ERA_precip_downcaled(:,count+1:count+4) = repmat(0,size(tile.PARA.latitude,1), 4);
                    
                end
                count = count + 4;
                
            end
            forcing.PARA.heatFlux_lb = 0.05;
            forcing.PARA.airT_height = 2;

        end
        
        
        function output = ncread_with_360_degree(forcing, fname, variable, time_index_start, time_number_of_values)
            
            if forcing.TEMP.max_target_lon < 360 %"normal" case, no need to change anything
                output = ncread(fname, variable, [ forcing.TEMP.lon_index_start forcing.TEMP.lat_index_start time_index_start], [forcing.TEMP.lon_index_end - forcing.TEMP.lon_index_start+1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 time_number_of_values], [1 1 1]);
                
            else  %this part is not yet tested, it only becomes effective at the 0 degree longitude -> the only affected PF seems to be in the Pyrenees, which is fairly irrelevant
                forcing.TEMP.lon_index_end=forcing.TEMP.lon_index_end-1; %read to lon = 359.75 degree
                output = ncread(fname, variable, [ forcing.TEMP.lon_index_start forcing.TEMP.lat_index_start time_index_start], [forcing.TEMP.lon_index_end - forcing.TEMP.lon_index_start+1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 time_number_of_values], [1 1 1]);
                output0degree = ncread(fname, variable, [ 1 forcing.TEMP.lat_index_start time_index_start], [1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 time_number_of_values], [1 1 1]); %read first "line" of longitude = 0 degree
                output = cat(1, output, output0degree); %concat in dimension1 (longitude)
            end
            
        end
        
        function output = ncread_with_360_degree_pl(forcing, fname, variable, time_index_start, time_number_of_values, pl)
            
            %pl = 1 1000; pl = 2 700; pl = 3 500; pl = 4 200
            
            if forcing.TEMP.max_target_lon < 360 %"normal" case, no need to change anything
                output = ncread(fname, variable, [ forcing.TEMP.lon_index_start forcing.TEMP.lat_index_start pl time_index_start], [forcing.TEMP.lon_index_end - forcing.TEMP.lon_index_start+1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 1 time_number_of_values], [1 1 1 1]);
            else  %this part is not yet tested, it only becomes effective at the 0 degree longitude -> the only affected PF seems to be in the Pyrenees, which is fairly irrelevant
                forcing.TEMP.lon_index_end=forcing.TEMP.lon_index_end-1; %read to lon = 359.75 degree
                output = ncread(fname, variable, [ forcing.TEMP.lon_index_start forcing.TEMP.lat_index_start pl time_index_start], [forcing.TEMP.lon_index_end - forcing.TEMP.lon_index_start+1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 1 time_number_of_values], [1 1 1 1]);
                output0degree = ncread(fname, variable, [ 1 forcing.TEMP.lat_index_start pl time_index_start], [1 forcing.TEMP.lat_index_end - forcing.TEMP.lat_index_start+1 1 time_number_of_values], [1 1 1 1]); %read first "line" of longitude
                output = cat(1, output, output0degree); %concat in dimension1 (longitude)
            end
        end
        
    end
end