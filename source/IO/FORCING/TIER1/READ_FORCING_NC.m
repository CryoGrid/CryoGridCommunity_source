%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class to read forcing data provided in NetCDF format
% 
% Robin B. Zweigel, October 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef READ_FORCING_NC < READ_FORCING_base
    
    properties
        
    end
    
    methods
% MAKE NEW SECTION OF READ FUNCTIONS FOR NEW DATA SOURCES
        
        
%---------------------- ERA5 read functions  ------------------------------
        function forcing = read_NC_ERA5_vars(forcing, path)
            variables = fields(forcing.DATA);
            
            for i = 1:length(variables)
                get_function = str2func(['ERA5_get_' variables{i}]);
                forcing.DATA.(variables{i}) = get_function(forcing, path);
            end
            
        end
        
        function [data, times] = read_NC_ERA5_vars_old(forcing, path, variables)
            
            for i=1:size(variables,1)
                data.(variables{i,1}) = squeeze(ncread([path variables{i,1} '.nc'], variables{i,1}));
            end
            
            temp.time = ncread([path variables{1,1} '.nc'], 'time');
            temp.info = ncinfo([path variables{1,1} '.nc']);
            Index1 = find(contains({temp.info.Variables.Name},'time'));
            Index2 = find(contains({temp.info.Variables(Index1).Attributes.Name},'units'));
            reference = split(temp.info.Variables(Index1).Attributes(Index2).Value);
            reference = reference(3);
            reference_date = split(reference,"-");
            reference_date = str2double(reference_date);
            times = datenum(reference_date(1),reference_date(2),reference_date(3)) + double(temp.time)./24;
            
        end
        
        function timeForcing = ERA5_get_timeForcing(forcing, path)
            temp.time = ncread([path 't2m.nc'], 'time');
            temp.info = ncinfo([path 't2m.nc']);
            Index1 = find(contains({temp.info.Variables.Name},'time'));
            Index2 = find(contains({temp.info.Variables(Index1).Attributes.Name},'units'));
            reference = split(temp.info.Variables(Index1).Attributes(Index2).Value);
            reference = reference(3);
            reference_date = split(reference,"-");
            reference_date = str2double(reference_date);
            timeForcing = datenum(reference_date(1),reference_date(2),reference_date(3)) + double(temp.time)./24;
        end
        
        function Tair = ERA5_get_Tair(forcing, path)
            t2m = squeeze(ncread([path 't2m.nc'], 't2m'));
            Tair = t2m - forcing.CONST.Tmfw;
        end
        
        function wind = ERA5_get_wind(forcing, path)
            u10 = squeeze(ncread([path 'u10.nc'], 'u10'));
            v10 = squeeze(ncread([path 'v10.nc'], 'v10'));
            wind = sqrt(u10.^2 + v10.^2);
        end
        
        function Sin = ERA5_get_Sin(forcing, path)
            ssrd = squeeze(ncread([path 'ssrd.nc'], 'ssrd'));
            Sin = [0; 0; ssrd./3600];
        end
        
        function Lin = ERA5_get_Lin(forcing, path)
            strd = squeeze(ncread([path 'strd.nc'], 'strd'));
            Lin = [0; 0; strd./3600];
        end
        
        function p = ERA5_get_p(forcing, path)
            sp = squeeze(ncread([path 'sp.nc'], 'sp'));
            p = sp;
            % Maybe make swich in case surface pressure is not provided????
            % RBZ oct. 2022
        end
        
        function q = ERA5_get_q(forcing, path)
            d2m = squeeze(ncread([path 'd2m.nc'], 'd2m'));
            q = 0.622.*(double(forcing.DATA.Tair<0).*satPresIce(forcing, d2m) + double(forcing.DATA.Tair>=0).*satPresWater(forcing, d2m))./ forcing.DATA.p;
        end
        
        function precip = ERA5_get_precip(forcing, path)
            tp = squeeze(ncread([path 'tp.nc'], 'tp'));
            precip = [0;0; tp];
        end
        
        function S_TOA = ERA5_get_S_TOA(forcing, path)
            tisr = squeeze(ncread([path 'tisr.nc'], 'tisr'));
            S_TOA = [0;0; tisr./3600];
        end
        
        
        %----------read functions for TopoScale------------------
        function forcing = load_ERA_slice(forcing)
            
            disp('load nc-files')
            
            if strcmp(forcing.PARA.time_resolution_input, 'quarter')
                quarter = forcing.TEMP.current_quarter;
            elseif strcmp(forcing.PARA.time_resolution_input, 'month')
                month = forcing.TEMP.current_month;
            end
            year = forcing.TEMP.current_year;
                        
            if strcmp(forcing.PARA.time_resolution_input, 'quarter')
                ncf=[forcing.PARA.forcing_path forcing.PARA.nc_folder 'surf_q' num2str(quarter) '_' num2str(year) '.nc'];
            elseif strcmp(forcing.PARA.time_resolution_input, 'month')
                ncf=[forcing.PARA.forcing_path forcing.PARA.nc_folder 'surf_m' num2str(month) '_' num2str(year) '.nc'];
            elseif strcmp(forcing.PARA.time_resolution_input, 'year')
                ncf=[forcing.PARA.forcing_path forcing.PARA.nc_folder 'surf_' num2str(year) '.nc'];                
            end
                        
            t = ncread(ncf,'time');
            target_timestep = (double(t(2))-double(t(1)))./24;
            
            if strcmp(forcing.PARA.time_resolution_input, 'quarter')
                era.t=datenum(year,(quarter-1).*3+1, 1):target_timestep:datenum(year, quarter.*3+1,1)-target_timestep; % Target time steps
            elseif strcmp(forcing.PARA.time_resolution_input, 'month')
                era.t=datenum(year,month, 1):target_timestep:datenum(year, month+1,1)-target_timestep; % Target time steps
            elseif strcmp(forcing.PARA.time_resolution_input, 'year')
                era.t=datenum(year,1, 1):target_timestep:datenum(year+1,1,1)-target_timestep; % Target time steps
            end

            %only needed during initialization
            era.lon=ncread(ncf,'longitude');
            era.lat=ncread(ncf,'latitude');
            
            % [era.P,era.SW,era.LW,era.ps,era.T2,era.Td2,era.u10,era.v10]=...
            %     deal(nan(numel(era.lon),numel(era.lat),numel(era.t),'single'));
            
            ncfz=[forcing.PARA.forcing_path forcing.PARA.nc_folder 'gp_surf.nc'];
            tmp=ncread(ncfz,'z');
            era.Zs=tmp(:,:,1)./9.81;

            u10=ncread(ncf,'u10');
            era.u10=u10;
            v10=ncread(ncf,'v10');
            era.v10=v10;
            Td2=ncread(ncf,'d2m');
            era.Td2=Td2;
            T2=ncread(ncf,'t2m');
            era.T2=T2;
            ps=ncread(ncf,'sp');
            era.ps=ps;
            SW=ncread(ncf,'ssrd');
            era.SW=SW./(60.^2);
            LW=ncread(ncf,'strd');
            era.LW=LW./(60.^2);
            P=ncread(ncf,'tp');
            era.P=P.*1e3; % month/hour to mm/hour
            
            p = forcing.PARA.top_pressure_level:25:forcing.PARA.bottom_pressure_level;%[800 825 850 875 900 925 950 975 1000]; era.p=p.*1e2; % Store in Pascals.
            era.p=[];
            
            [era.Z,era.q,era.T,era.u,era.v]=...
                deal(nan(numel(era.lon),numel(era.lat),1,numel(era.t),'single'));
            
            % Loop over pressure levels
            ind = 1;
            for l=1:numel(p)
                if strcmp(forcing.PARA.time_resolution_input, 'quarter')
                    ncf=[forcing.PARA.forcing_path forcing.PARA.nc_folder 'plev_q' num2str(quarter) '_p' num2str(p(l)) '_' num2str(year) '.nc'];
                elseif strcmp(forcing.PARA.time_resolution_input, 'month')
                    ncf=[forcing.PARA.forcing_path forcing.PARA.nc_folder 'plev_m' num2str(month) '_p' num2str(p(l)) '_' num2str(year) '.nc'];
                elseif strcmp(forcing.PARA.time_resolution_input, 'year')
                    ncf=[forcing.PARA.forcing_path forcing.PARA.nc_folder 'plev_p' num2str(p(l)) '_' num2str(year) '.nc'];
                end
                if exist(ncf)==2
                    Z=ncread(ncf,'z')./9.81;
                    era.Z(:,:,ind,:)=Z;
                    qp=ncread(ncf,'q'); % Named quarter counter "q" as well
                    era.q(:,:,ind,:)=qp;
                    T=ncread(ncf,'t');
                    era.T(:,:,ind,:)=T;
                    u=ncread(ncf,'u');
                    era.u(:,:,ind,:)=u;
                    v=ncread(ncf,'v');
                    era.v(:,:,ind,:)=v;
                    era.p=[era.p p(l).*100];
                    ind = ind+1;
                end
            end

            % Scale output.
            era.wind_sf=1e-2;
            era.q_sf=1e-6;
            era.ps_sf=1e2;
            era.rad_sf=1e-1;
            era.T_sf=1e-2;
            era.P_sf=1e-2;
            
            era.u=int16(era.u./era.wind_sf);
            era.u10=int16(era.u10./era.wind_sf);
            era.v10=int16(era.v10./era.wind_sf);
            era.v=int16(era.v./era.wind_sf);
            era.q=uint16(era.q./era.q_sf);
            era.ps=uint16(era.ps./era.ps_sf);
            era.SW=uint16(era.SW./era.rad_sf);
            era.LW=uint16(era.LW./era.rad_sf); % Do this in the time loop.
            era.T=int16((era.T-273.15)./era.T_sf);
            era.T2=int16((era.T2-273.15)./era.T_sf);
            era.Td2=int16((era.Td2-273.15)./era.T_sf);
            era.P=uint16(era.P./era.P_sf);
            era.Z=int16(era.Z);
            
            forcing.TEMP.era = era;
        end
        
        function forcing = load_ERA_sl_slice(forcing)
            
            disp('load nc-files')
            
            if strcmp(forcing.PARA.time_resolution_input, 'quarter')
                quarter = forcing.TEMP.current_quarter;
            elseif strcmp(forcing.PARA.time_resolution_input, 'month')
                month = forcing.TEMP.current_month;
            end
            year = forcing.TEMP.current_year;
                        
            if strcmp(forcing.PARA.time_resolution_input, 'quarter')
                ncf=[forcing.PARA.forcing_path forcing.PARA.nc_folder 'surf_q' num2str(quarter) '_' num2str(year) '.nc'];
            elseif strcmp(forcing.PARA.time_resolution_input, 'month')
                ncf=[forcing.PARA.forcing_path forcing.PARA.nc_folder 'surf_m' num2str(month) '_' num2str(year) '.nc'];
            elseif strcmp(forcing.PARA.time_resolution_input, 'year')
                ncf=[forcing.PARA.forcing_path forcing.PARA.nc_folder 'surf_' num2str(year) '.nc'];                
            end
                        
            t = ncread(ncf,'time');
            target_timestep = (double(t(2))-double(t(1)))./24;
            
            if strcmp(forcing.PARA.time_resolution_input, 'quarter')
                era.t=datenum(year,(quarter-1).*3+1, 1):target_timestep:datenum(year, quarter.*3+1,1)-target_timestep; % Target time steps
            elseif strcmp(forcing.PARA.time_resolution_input, 'month')
                era.t=datenum(year,month, 1):target_timestep:datenum(year, month+1,1)-target_timestep; % Target time steps
            elseif strcmp(forcing.PARA.time_resolution_input, 'year')
                era.t=datenum(year,1, 1):target_timestep:datenum(year+1,1,1)-target_timestep; % Target time steps
            end

            %only needed during initialization
            era.lon=ncread(ncf,'longitude');
            era.lat=ncread(ncf,'latitude');
            
            % [era.P,era.SW,era.LW,era.ps,era.T2,era.Td2,era.u10,era.v10]=...
            %     deal(nan(numel(era.lon),numel(era.lat),numel(era.t),'single'));
            
            ncfz=[forcing.PARA.forcing_path forcing.PARA.nc_folder 'gp_surf.nc'];
            tmp=ncread(ncfz,'z');
            era.Zs=tmp(:,:,1)./9.81;

            u10=ncread(ncf,'u10');
            era.u10=u10;
            v10=ncread(ncf,'v10');
            era.v10=v10;
            Td2=ncread(ncf,'d2m');
            era.Td2=Td2;
            T2=ncread(ncf,'t2m');
            era.T2=T2;
            ps=ncread(ncf,'sp');
            era.ps=ps;
            SW=ncread(ncf,'ssrd');
            era.SW=SW./(60.^2);
            LW=ncread(ncf,'strd');
            era.LW=LW./(60.^2);
            P=ncread(ncf,'tp');
            era.P=P.*1e3; % m/hour to mm/hour

            % Scale output.
            era.wind_sf=1e-2;
            era.q_sf=1e-6;
            era.ps_sf=1e2;
            era.rad_sf=1e-1;
            era.T_sf=1e-2;
            era.P_sf=1e-2;
            
            era.u10=int16(era.u10./era.wind_sf);
            era.v10=int16(era.v10./era.wind_sf);
            era.ps=uint16(era.ps./era.ps_sf);
            era.SW=uint16(era.SW./era.rad_sf);
            era.LW=uint16(era.LW./era.rad_sf); % Do this in the time loop.
            era.T2=int16((era.T2-273.15)./era.T_sf);
            era.Td2=int16((era.Td2-273.15)./era.T_sf);
            era.P=uint16(era.P./era.P_sf);
%             era.Z=int16(era.Z);
            
            forcing.TEMP.era = era;
        end
    end
    
    
end
