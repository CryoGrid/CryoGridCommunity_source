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
                %%NC added 
                if strcmp(section{i,1}, 'leaching_rate')
                    forcing.PARA.leaching_rate = section{i,2};
                end
                
                if strcmp(section{i,1}, 'total_tracer')
                    forcing.PARA.total_tracer = section{i,2};
                end
                
                if strcmp(section{i,1}, 'infiltration_cutoff')
                    forcing.PARA.infiltration_cutoff = section{i,2};
                end
                
                if strcmp(section{i,1}, 'infiltration')
                    forcing.PARA.module.infiltration = section{i,2};
                end
                
                if strcmp(section{i,1}, 'xice')
                    forcing.PARA.module.xice = section{i,2};
                end
                
                if strcmp(section{i,1}, 'lateral')
                    forcing.PARA.module.lateral = section{i,2};
                end
                
                if strcmp(section{i,1}, 'exchange_heat')
                    forcing.PARA.module.exchange_heat = section{i,2};
                end 
                
                if strcmp(section{i,1}, 'exchange_water')
                    forcing.PARA.module.exchange_water = section{i,2};
                end
                
                if strcmp(section{i,1}, 'exchange_snow')
                    forcing.PARA.module.exchange_snow = section{i,2};
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