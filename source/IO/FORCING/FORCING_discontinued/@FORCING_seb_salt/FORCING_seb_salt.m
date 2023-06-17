% forcing data


classdef FORCING_seb_salt
    properties
        DATA  %all data
        TEMP  %at each timestep
        PARA
        STATUS %forcing data suitable for the modules that are to be run -> can be used
    end


    methods

        function forcing = set_parameters(forcing, filename, rain_fraction, snow_fraction, altitude, salt_concentration) %salt_c can be given as parameter, but gets overwritten when time series is available in the focing mat-file
            forcing.PARA.filename = filename;
            forcing.PARA.rain_fraction = rain_fraction;
            forcing.PARA.snow_fraction = snow_fraction;
            forcing.PARA.altitude = altitude;
            forcing.PARA.saltConcForcing = salt_concentration;
        end

       % function forcing = load_forcing_from_mat(forcing);
       % function forcing = interpolate_forcing(t, forcing);

    end
end
