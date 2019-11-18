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
           

           
           %calculate time in days for the main programm
           forcing.PARA.start_time = forcing.PARA.startForcing .* 365.25;
           forcing.PARA.end_time = forcing.PARA.endForcing .*365.25;

       end

       function forcing = interpolate_forcing(t, forcing) %t comes in days!
           t = t/365.25; %convert to t in years
           forcing.TEMP.TForcing = interp1(forcing.DATA.timeForcing, forcing.DATA.TForcing, t);
           %forcing.TEMP.surfaceState = interp1(forcing.DATA.timeForcing, forcing.DATA.surfaceState, t, 'nearest');

           forcing.TEMP.saltConcForcing = interp1(forcing.DATA.timeForcing, forcing.DATA.saltConcForcing, t, 'nearest');

           %> For the forcing of the salt, we use the surfaceState
           %> (subaerial, subglacial, submarine) as a flag for the salt flux
           %> To avoid rapid changes, we interpolate lineraly if the change
           %> is from or to a submarine phase.
%            stateBefore = interp1(forcing.DATA.timeForcing, forcing.DATA.surfaceState, t, 'previous');
%            stateAfter = interp1(forcing.DATA.timeForcing, forcing.DATA.surfaceState, t, 'next');
%            if stateBefore == 0 || stateAfter == 0 %if one of the states is submarine
%                 forcing.TEMP.surfaceState = interp1(forcing.DATA.timeForcing, forcing.DATA.surfaceState, t);
%            else
%                 forcing.TEMP.surfaceState = interp1(forcing.DATA.timeForcing, forcing.DATA.surfaceState, t, 'nearest');
%            end

           %This is from the old getDerivative/testMexSubseaPF
%            	%What is happening here? Is this relevant for the upper boundary?
%             %Should this happen in ground.get_boundary_condition_u?
%             deltaT = TForcing(2) - TForcing(1);
%             factor = floor((t/day_sec - TForcing(1))/deltaT);
%             T_u = TForcing(factor + length(TForcing)) + ...
%                 (Tsurf(factor + 1 + length(TForcing)) - TForcing(factor + length(TForcing))) * ...
%                 (t/day_sec - TForcing(factor))/deltaT;
       end

    end
end
