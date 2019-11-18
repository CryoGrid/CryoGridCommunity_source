% forcing data


classdef FORCING_seb
    properties
        DATA  %all data
        TEMP  %at each timestep
        PARA
        STATUS %forcing data suitable for the modules that are to be run -> can be used 
    end
    
    
    methods
        
        function forcing = initalize_from_file(forcing, section)
            for i=1:size(section,1)
                if strcmp(section{i,1}, 'filename')
                    forcing.PARA.filename = section{i,2};
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