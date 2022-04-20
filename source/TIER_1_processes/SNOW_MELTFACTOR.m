%========================================================================
% CryoGrid TIER1 library class for functions related to snow cover
% melting using simple melt factor.
%
% The melt factor is in the range of observed values (2-12 mmSWE/day) and 
% scaled by the length of day and maximum sun angle (derived from the date
% and the latitude).
% The melt rate is then given by the product of the melt factor and the air
% temperature in given time step.
% Knowing the melt rate, the energy needed to obtain this melt rate within 
% the current timestep is calculated based on the temperature and density 
% of the uppermost snow cell, and timestep is chosen such that within the
% timestep no more than half of a cell is melted.
%
% Citations:
% Hock, R. (2003): Temperature index melt modelling in mountain areas.
%     Journal of Hydrology 282, pp 104-115.
% Obu et al. (2019): Northern Hemisphere permafrost map based on TTOP
%     modelling for 2000-2016 at 1 km2 scale. Earth-Science Reviews 193,
%     pp. 299-316.
%
% Implemented by:
% T. Ingeman-Nielsen, S. Westermann, December 2020
%========================================================================

classdef SNOW_MELTFACTOR < BASE
    
    methods
        
        %--boundary conditions--------
        
        function snow = get_boundary_condition_SNOW_meltFactor(snow, tile)
            forcing = tile.FORCING;
            % Calculates the meltfactor and assoicated energy contribution 
            % We assume a melt factor (in range 2-12 mm/day) which is
            % scaled by the length of day and maximum sun angle.
            %
            % We assume melt rate is in snow water equivalents, that
            % all snow melts from the top of the snowpack, and that never
            % more than one cell melts in a particular time step.
            
            L_f = snow.CONST.L_f;   % volumetric latent heat of fusion, freezing [J/m3]
            c_i = snow.CONST.c_i;  % volumetric heat capacity of ice [J/m3/K]
            
            % calculate the melt factor [mm/day/degC] from timestamp and location
            snow.TEMP.melt_factor = SNOW_MELTFACTOR.MeltFactorsFromDayLenAndSunAngle(tile.t, tile.PARA.latitude);
                        
            % calculate the melt rate [m3/s] for the current time interval
            % melt_rate is given in snow water equivalents
            snow.TEMP.melt_rate = (snow.TEMP.melt_factor./1000./snow.CONST.day_sec) .* ...
                (forcing.TEMP.Tair - snow.PARA.melt_threshold) .* snow.STATVAR.area(1,1);
            snow.TEMP.melt_rate = max(0,snow.TEMP.melt_rate);
            % 1000 mm/m is the length conversion factor
            % snow.CONST.day_sec is the factor used to convert from days to seconds
            % snow.PARA.melt_threshold is the threshold air temperature above which snow melt occurs [degC]
            
            % calculate the melt energy rate
            snow.TEMP.melt_energy = snow.TEMP.melt_rate .* (-1*min(0, snow.STATVAR.T(1))*c_i + L_f); % [J/s] 
            % only calculated for top cell, since snow melt happens from
            % top of snow pack, and the timestep is chosen so that never
            % more than one cell is melted.
            
        end
        
        
        %--timesteps------------------
        
        % Time step modification is handled by the standard methods of the
        % TIER1 SNOW class.
        
        %--diagnostic step------------
        
        
        %--triggers-------------------
        
        
        %--other methods--------------
        
    end
    
    methods(Static)
        function [mfs] = MeltFactorsFromDayLenAndSunAngle(timestamp, latitude)
            % Function calculates average melt factors (MFs) for a specific day
            % based on sun culmination angle and day length. For this, it uses "day_length" and
            % "solarCulminationAngle" functions. The MFs are linearly scaled between a
            % product of culmination angle and day length.
            % 
            % Citation:
            % XXXXXXXXXXX

            % maximum and minimum MFs to which MFs are scaled 
            mfmin = 2;
            mfmax = 12;

            % maximum product between culmination angle and day length/2 in degrees.
            maxProduct = 10169.237;

            % solar culmination angle
            curcul = SNOW_MELTFACTOR.solarCulminationAngle(timestamp, latitude);
            % only positive
            curcul(curcul < 0) = 0;

            % day of the year according to winter solstice
            dlDoy = SNOW_MELTFACTOR.doy_int(timestamp)+11;
            % day lengths in hours
            curlen = SNOW_MELTFACTOR.day_length(dlDoy, latitude);
            curlen(curlen < 0) = 0;
            % scaling hours to degrees. 
            curlen_deg = (curlen/24)*360;

            % multiplying culmination angle with half of the daylength degrees. The
            % product equals approx. area of sun path above horizont.
            daylculm = curcul.*(curlen_deg/2);
            % converting to fraction based on maximum šroduct
            daylculm_frac = daylculm/maxProduct;

            % scaling MFs to maximum and minimum DDF
            mfs = daylculm_frac*(mfmax-mfmin)+mfmin;
        end
        
        
        function [hours, b] = day_length(Day, Latitude)
            % [hours, b] = day_length(Day, Latitude)
            % This calculates the number of hours (hours) and fraction of the day (b) in %daylight.
            % 
            % Inputs:
            % Day:       day of the year, counted starting with the day of 
            %                the December solstice in the first year of a 
            %                Great Year.
            % Latitude:  latitude in degrees, North is positive, 
            %                South is negative
            %
            % Calculations are per Herbert Glarner's formulae which do not take into account refraction, 
            % twilight, size of the sun, etc. 
            % (http://herbert.gandraxa.com/herbert/lod.asp but be careful about inconsistencies in radians/degrees).
            %
            % Copyright (c) 2015, Travis Wiens
            % All rights reserved.

            if nargin < 1
              Day=0;   % default to first day of year
            end
            if nargin < 2
              Latitude = (-(36+51/60)); % default to 36d51'
            end

            Axis = 23.439*pi/180;

            j = pi/182.625; % constant (radians)
            m = 1-tan(Latitude*pi/180).*tan(Axis*cos(j*Day));

            m(m>2) = 2; % saturate value for artic
            m(m<0) = 0;

            b = acos(1-m)/pi; % fraction of the day the sun is up

            hours = b*24; % hours of sunlight
        end
        

        function doy = doy_int(date)
            % create integers of day of year (1-366), from datenum

            f = zeros(size(date));
            [y,m,d,~,~,~] = datevec(date);
            doy = datenum(y,m,d,f,f,f) - datenum(y,f,f,f,f,f);
        end

        
        function [culmination_angle, declination_angle] = solarCulminationAngle(timestamp, latitude)
            % calculates the solar culmination and declination angle
            %
            % The culmination angle is the angle of the sun above the 
            % observers local horizon (assuming a spherical earth).
            % The declination angle is the angle of the sun above the
            % equatorial plane of the earth.
            
            %day of the year
            n = SNOW_MELTFACTOR.doy_int(timestamp);

            %declination angle formula
            declination_angle = 23.45 * sind(360/365*(284+n));
            culmination_angle = 90-latitude+declination_angle;
        end

        
    end
end

