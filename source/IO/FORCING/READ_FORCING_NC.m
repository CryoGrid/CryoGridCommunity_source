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
            wind = sqrt(u10.^2 + u10.^2);
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
            q = (double(forcing.DATA.Tair<0).*satPresIce(forcing, d2m) + double(forcing.DATA.Tair>=0).*satPresWater(forcing, d2m))./ forcing.DATA.p;
        end
        
        function precip = ERA5_get_precip(forcing, path)
            tp = squeeze(ncread([path 'tp.nc'], 'tp'));
            precip = [0;0; tp];
        end
        
        function S_TOA = ERA5_get_S_TOA(forcing, path)
            tisr = squeeze(ncread([path 'tisr.nc'], 'tisr'));
            S_TOA = [0;0; tisr./3600];
        end
        
         
    end
    
    
end
