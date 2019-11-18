function forcing = generateForcing_fromData(forcing)
%   generateForcing accepts (lat, lon, zsb) for the location and
%   returns the timestamp, the duration
%   of inundation, surface temperature and an index defining the surface
%   (submarine, subaerial, subglacial).

lat = forcing.PARA.latitude;
lon = forcing.PARA.longitude;
zsb = forcing.PARA.altitude;

timeForcing = [forcing.PARA.startForcing:forcing.PARA.dtForcing:forcing.PARA.endForcing]';
forcing.DATA.timeForcing = timeForcing;
saltConcForcing = timeForcing.*0; %"outer" salt concentration is zero for subaerial and subglacial conditions
coverage = timeForcing.*0;    % positive values are thickness of ice-coverage
                            % negative values are depth of covering sea
                            % (over current elevation)
                            % 0 accounts for no coverage

%%  define sea level history
switch forcing.PARA.SL_no
    case 1% Lambeck et al.
        load('forcing/subseaInput.mat','elevation','nominal_elevation','time')
        seaLevel = interp1([time(1:325,1)*1000; 50000], [nominal_elevation(1:325,1); nominal_elevation(325,1)], -timeForcing);
    case 2 % Walbroeck et al.
        load('forcing/sealevel_Waelbroeck.mat')  % starting at -429.5 kyrs BP
        time_W = [0 ; time_W]; RSL_W = [0 ; RSL_W]; % add RSL=0 for t=0
        seaLevel = interp1(time_W*1000,RSL_W,timeForcing);
    case 3 % Grant et al.
        load('forcing/sealevel_Grant_prob.mat')  % starting at -492 kyrs BP
        seaLevel = interp1(time_G*1000,RSL_G-RSL_G(1),timeForcing); % set RSL=0 at year t=0 (i.e. substract -RSL_G(1))

    case 666 %testing without inputdata: linear rise from -20 to 0 over time
        seaLevel = linspace(-20,0,length(timeForcing));

    otherwise
        disp('Invalid switch for sea level history')
end
forcing.DATA.seaLevel = seaLevel;


%%  define air temperature forcing
switch forcing.PARA.TF_no
    case 1 % use GISP-2 data
        load('forcing/GISP2.mat', 'GISP2_T')
        baseT = -5;   % air temperature offset for GISP data
        for i=size(GISP2_T,1):-1:2
            if GISP2_T(i,1)==GISP2_T(i-1,1)
                GISP2_T(i,:)=[];
            end
        end
        deltaT=interp1([0; GISP2_T(:,1)*1000; 50000], [GISP2_T(1,2); GISP2_T(:,2); GISP2_T(end,2)], -timeForcing);
        deltaT=deltaT-deltaT(end,1);
        T = baseT+deltaT;
        forcingData = T;
    case 2 % use CLIMBER-2 data
%         if(EDF.GR==1)
%             load forcing/CLIMBER2/Data_indices_lat_lon_C2_all_25km  % get index data (lat&lon) for CLIMBER-2 SAT grid
%         elseif(EDF.GR==2)
            load forcing/CLIMBER2/Data_indices_lat_lon_C2_all_12p5km  %
%         elseif(EDF.GR==3)
%             load InputData/CLIMBER2/Data_indices_lat_lon_C2_all_6p25km  %
%         end
        % load InputData/CLIMBER2/SAT_ANN_C2ip_500k         % get CLIMBER-2 SAT matrix (interpolated to grid 1.5?lon [-178.9 179.6], 0.75?lat [55N 85N], annually averaged, 100 yrs timeslices)
        load forcing/CLIMBER2/SATcor_ANN_C2ip_500k          % CLIMBER-2 SAT, corrected for thermal offset
        %t1 = length(time_C2ip)+(forcing.PARA.startForcing/forcing.PARA.dtForcing);
        %t2 = length(time_C2ip); % get index of SAT time series
        this_ind_lat = find(abs(lat_C2ip-lat)==min(abs(lat_C2ip-lat)), 1, 'first');
        this_ind_lon = find(abs(lon_C2ip2-lon)==min(abs(lon_C2ip2-lon)), 1, 'first');
        %SAT_C2 = squeeze(SAT_ANNipt2(this_ind_lon,this_ind_lat,t1:t2)); % time series of CLIMBER-2 SAT for chosen grid cell
        SAT_C2 = squeeze(SAT_ANNipt2(this_ind_lon,this_ind_lat,:)); % time series of CLIMBER-2 SAT for chosen grid cell
        SAT_C2_interp = interp1(time_C2ip, SAT_C2, timeForcing); %interpolate linear to given time steps
        clear SAT_ANNipt2  % free memory
        forcingData = SAT_C2_interp;

    case 666 %testing without inputdata: linear rise from -10 to -1 over time
        forcingData = linspace(-10,-1,length(timeForcing));
    otherwise
        disp('Invalid switch TF for specifying temperature forcing')
end
forcing.DATA.airTemp = forcingData;

%% set temperature data for inundated sites to sea water temperature
ind_inundation = find(seaLevel > zsb); % index of years of inundation
time_inundation = forcing.PARA.dtForcing * length(ind_inundation); % duration of inundation

for i=1:length(forcingData) % adapt forcing temperature for inundated sites
    if(zsb < seaLevel(i)) % site is inundated
        waterDepth = seaLevel(i) - zsb;    % depth water column
        if(waterDepth > 30) % below 30m T sea bottom equals T_freeze)
            T_seaWater = forcing.PARA.T_freeze;
        elseif(waterDepth <= 30 && waterDepth > 2) % linear scaling between 30m t0 2m water depth
            T_seaWater = 1./14 * (forcing.PARA.T_freeze/2*waterDepth - forcing.PARA.T_freeze);
        else  % between 2m and 0m T sea bottom equals 0�C
            T_seaWater = 0;
        end
        forcingData(i) = T_seaWater;

        coverage(i) = -waterDepth;

        saltConcForcing(i) = forcing.PARA.benthicSalt;
    end
end

% set temperature data for (warm-based) ice sheet covered sites to bottom ice sheet temperature
if(forcing.PARA.IS > 0) % use of CLIMBER-2 ice sheet data
    load forcing/CLIMBER2/Ice_Thickness_C2ipt_500kyrs   % get CLIMBER-2 ice sheet thickness matrix (same grid resolution as SAT matrix
    if ~exist('this_ind_lat', 'var')
        this_ind_lat = find(abs(lat_C2ip-lat)==min(abs(lat_C2ip-lat)), 1, 'first');
        this_ind_lon = find(abs(lon_C2ip2-lon)==min(abs(lon_C2ip2-lon)), 1, 'first');
    end
    %ITH = squeeze(ITHipt(this_ind_lon,this_ind_lat,t1:t2));
    ITH = squeeze(ITHipt(this_ind_lon,this_ind_lat,:)); % time series of CLIMBER-2 ITH for chosen grid cell
    ITH = interp1(time_C2ip, ITH, timeForcing);
    clear ITHipt  % free memory
    ind_IceCover = find(ITH > forcing.PARA.IS);
    forcingData(ind_IceCover) = forcing.PARA.T_IceSheet; % set temperature to ice sheet bottom temperatures for sites with a ice sheet coverage larger than S_IS
    ind_IS = ITH > forcing.PARA.IS; % index of years of IceSheet coverage thicker than S_IS

    coverage(ind_IceCover) = ITH(ind_IceCover);
end


surfaceState=timeForcing.*0 + 1;      % define surface state as 1 everywhen (subaerial==1)
surfaceState(ind_inundation) = 0;   % set inundation periods to 0 (submarine==0)
if exist('ind_IS', 'var')
    surfaceState(ind_IS) = -1;         % set subglacial periods to -1 (subglacial ==-1)
end


forcing.DATA.TForcing = forcingData;
forcing.DATA.coverage = coverage;
forcing.DATA.time_inundation = time_inundation;
forcing.DATA.surfaceState = surfaceState;
if forcing.PARA.saltForcingSwitch == 0
    forcing.DATA.saltConcForcing = zeros(size(saltConcForcing));
else
    forcing.DATA.saltConcForcing = saltConcForcing;
end
