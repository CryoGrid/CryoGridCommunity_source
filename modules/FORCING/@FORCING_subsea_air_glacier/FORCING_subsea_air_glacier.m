% forcing data


classdef FORCING_subsea_air_glacier
    properties
        DATA  %all data
        TEMP  %at each timestep
        PARA
        STATUS %forcing data suitable for the modules that are to be run -> can be used
    end


    methods
        function forcing = initalize_from_file(forcing, section)
           forcing.PARA.SL_no = [];  
           forcing.PARA.TF_no = [];
           forcing.PARA.T_freeze = [];
           forcing.PARA.T_IceSheet = [];
           forcing.PARA.IS = [];
           forcing.PARA.benthicSalt = [];
           forcing.PARA.saltForcingSwitch = [];
           forcing.PARA.latitude = [];
           forcing.PARA.longitude = [];
           forcing.PARA.altitude = [];
           forcing.PARA.heatFlux_lb = [];
           forcing.PARA.startForcing = [];
           forcing.PARA.dtForcing = [];
           forcing.PARA.endForcing = [];
            
           forcing.PARA.domain_depth = 10000;
             
           variables = fieldnames(forcing.PARA);
           for i=1:size(variables,1)
               for j=1:size(section,1)
                  if strcmp(variables{i,1}, section{j,1})
                      forcing.PARA.(variables{i,1}) = section{j,2};
                  end
               end
           end

        end
          
       function forcing = load_forcing_from_mat(forcing)
           addpath(genpath('forcing'))
           
           %overwrite given altitude and heat flux by values from EaseGrid
           %and Davies et. al. respectively
           forcing.PARA.altitude = getElevation(forcing.PARA.latitude, forcing.PARA.longitude);
           forcing.PARA.heatFlux_lb = getQ_Davies(forcing.PARA.latitude, forcing.PARA.longitude);
           
           if forcing.PARA.SL_no == 666
               forcing = generateForcing_testing(forcing);
           else
               forcing = generateForcing_fromData(forcing);
           end
           
            %initialize elevation with modern-day altitude
            forcing.TEMP.elev = forcing.PARA.altitude;
           
           %calculate time in days for the main programm
           forcing.PARA.start_time = forcing.PARA.startForcing .* 365.25;
           forcing.PARA.end_time = forcing.PARA.endForcing .*365.25;

       end

       function forcing = interpolate_forcing(t, forcing) %t comes in days!
            t = t/365.25; %convert to t in years

            %get current sea level
            seaLevel = interp1(forcing.DATA.timeForcing, forcing.DATA.seaLevel, t);
            
            %get current altitude / upperPos
            elevation = forcing.TEMP.elev;
            
            %get glacial cover
            glacialCover = interp1(forcing.DATA.timeForcing, forcing.DATA.glacialCover, t);
            
            %get current forcing temperature depending on the surface state
            %Site is under water
            if elevation < seaLevel  % site is inundated
                waterDepth = seaLevel - elevation;    % depth water column
                if(waterDepth > 30) % below 30m T sea bottom equals T_freeze)
                    TForcing = forcing.PARA.T_freeze;
                elseif(waterDepth <= 30 && waterDepth > 2) % linear scaling between 30m t0 2m water depth
                    TForcing = 1./14 * (forcing.PARA.T_freeze/2*waterDepth - forcing.PARA.T_freeze);
                else  % between 2m and 0m T sea bottom equals 0�C
                    TForcing = 0;
                end
                saltConcForcing = forcing.PARA.benthicSalt;
                surfaceState = 0;
                
            elseif glacialCover > forcing.PARA.IS %if glacial cover is greater than treshold
                TForcing = forcing.PARA.T_IceSheet;
                saltConcForcing = 0;
                surfaceState = -1;
                
            else %subaerial conditions
                %get current air temperature
                TForcing = interp1(forcing.DATA.timeForcing, forcing.DATA.airTemp, t);
                saltConcForcing = 0;
                surfaceState = 1;
            end
    
            %update forcing struct
           forcing.TEMP.TForcing = TForcing;
           forcing.TEMP.saltConcForcing = saltConcForcing;
           forcing.TEMP.surfaceState = surfaceState;
       end

    end
end
