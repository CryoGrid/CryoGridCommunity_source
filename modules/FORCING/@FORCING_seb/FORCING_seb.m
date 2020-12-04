% forcing data


classdef FORCING_seb
    properties
        DATA  %all data
        TEMP  %at each timestep
        PARA
        STATUS %forcing data suitable for the modules that are to be run -> can be used 
    end
    
    
    methods
        
        function xls_out = write_excel(forcing)
            xls_out = {'FORCING','index',NaN,NaN;'FORCING_seb',1,NaN,NaN;NaN,NaN,NaN,NaN;'filename',NaN,NaN,NaN;'start_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the first timestamp of the forcing data set will be used';'end_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the last timestamp of the forcing data set will be used';'rain_fraction',1,'[-]','rainfall in forcing file multiplied by this number';'snow_fraction',1,'[-]','snowfall in forcing file multiplied by this number';'latitude',NaN,'[degree]','geographical coordinates';'longitude',NaN,'[degree]',NaN;'altitude',NaN,'[m]','a.s.l.';'domain_depth',100,'[m]','should match a GRID point, model domain extends to this depth';'heatFlux_lb',0.0500000000000000,'[W/m2]','geothermal heat flux';'airT_height',2,'[m]','height of air temperature';'FORCING_END',NaN,NaN,NaN};
        end
        
        function forcing = initalize_from_file(forcing, section)
            for i=1:size(section,1)
                if strcmp(section{i,1}, 'filename')
%                     forcing.PARA.filename = section{i,2};
%                     forcing.PARA.filename = 'Chukotka_tcc_T0_ERAint_1979_2019_trialSimone';
%                     forcing.PARA.filename = 'Nyurba_tcc_T0_ERAint_1979_2019_trialSimone';
%                      forcing.PARA.filename = 'Spasskaya_tcc_T0_ERAint_1979_2019_trialSimone';
                   forcing.PARA.filename = 'Mongolei_tcc_T0_ERAint_1979_2019_trialSimone';
                end
                
                
                if strcmp(section{i,1}, 'rain_fraction')
                    forcing.PARA.rain_fraction = section{i,2};
                end
                if strcmp(section{i,1}, 'snow_fraction')
                    forcing.PARA.snow_fraction = section{i,2};
                end
                if strcmp(section{i,1}, 'latitude')
                    forcing.PARA.latitude = section{i,2};
                end
                if strcmp(section{i,1}, 'longitude')
                    forcing.PARA.longitude = section{i,2};
                end
                if strcmp(section{i,1}, 'altitude')
                    forcing.PARA.altitude = section{i,2};
                end
                if strcmp(section{i,1}, 'domain_depth')
                    forcing.PARA.domain_depth = section{i,2};
                end
                if strcmp(section{i,1}, 'heatFlux_lb')
                    forcing.PARA.heatFlux_lb = section{i,2};
                end
                if strcmp(section{i,1}, 'airT_height')
                    forcing.PARA.airT_height = section{i,2};
                end
                if strcmp(section{i,1}, 'start_time')
                    forcing.PARA.start_time = section{i,2};
%                     forcing.PARA.start_time = '10.08.2015'; %start time for Batagai, ERA 5
                end
                if strcmp(section{i,1}, 'end_time')
                    forcing.PARA.end_time = section{i,2};
                end
                
                
            end
        end
        
        
        function forcing = set_parameters(forcing, filename, rain_fraction, snow_fraction, altitude)
            forcing.PARA.filename = filename;
            forcing.PARA.rain_fraction = rain_fraction;
            forcing.PARA.snow_fraction = snow_fraction;
            forcing.PARA.altitude = altitude;
        end
        
    end
end