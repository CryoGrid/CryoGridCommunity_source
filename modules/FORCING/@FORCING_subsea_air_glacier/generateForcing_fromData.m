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
        airTemp = T;
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
        airTemp = SAT_C2_interp;
    otherwise
        disp('Invalid switch TF for specifying temperature forcing')
end


%% define ice sheet coverage
if(forcing.PARA.IS > 0) % use of CLIMBER-2 ice sheet data
    load forcing/CLIMBER2/Ice_Thickness_C2ipt_500kyrs   % get CLIMBER-2 ice sheet thickness matrix (same grid resolution as SAT matrix
    if ~exist('this_ind_lat', 'var')
        this_ind_lat = find(abs(lat_C2ip-lat)==min(abs(lat_C2ip-lat)), 1, 'first');
        this_ind_lon = find(abs(lon_C2ip2-lon)==min(abs(lon_C2ip2-lon)), 1, 'first');
    end
    glacialCover = squeeze(ITHipt(this_ind_lon,this_ind_lat,:)); % time series of CLIMBER-2 ITH for chosen grid cell
    glacialCover = interp1(time_C2ip, glacialCover, timeForcing);
    clear ITHipt  % free memory
end




forcing.DATA.airTemp = airTemp;
forcing.DATA.seaLevel = seaLevel;
forcing.DATA.glacialCover = glacialCover;
