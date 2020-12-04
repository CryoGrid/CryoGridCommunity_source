function forcing = load_forcing_from_mat(forcing)

temp=load(['forcing/' forcing.PARA.filename], 'FORCING');

forcing.DATA.rainfall=temp.FORCING.data.rainfall.*forcing.PARA.rain_fraction;
forcing.DATA.snowfall=temp.FORCING.data.snowfall.*forcing.PARA.snow_fraction;
forcing.DATA.Tair = temp.FORCING.data.Tair;
forcing.DATA.Lin = temp.FORCING.data.Lin;
forcing.DATA.Sin = temp.FORCING.data.Sin;

forcing.DATA.Sin_dif = temp.FORCING.data.Sin_Irradiance_diffuse;
forcing.DATA.Sin_dir = temp.FORCING.data.Sin_Irradiance_direct;

forcing.DATA.q = temp.FORCING.data.q;
forcing.DATA.wind = temp.FORCING.data.wind;
forcing.DATA.timeForcing = temp.FORCING.data.t_span;
forcing.DATA.p = 101325; %temp.FORCING.data.pressure; % simone!

% Update spatial data if included in forcing file
if isfield(temp.FORCING.data,'z')
    forcing.PARA.altitude = round(temp.FORCING.data.z);
end
if isfield(temp.FORCING.data,'lon') &&isfield(temp.FORCING.data,'lat')
     forcing.PARA.longitude = temp.FORCING.data.lon;
     forcing.PARA.latitude  = temp.FORCING.data.lat;
end



if std(forcing.DATA.timeForcing(2:end,1)-forcing.DATA.timeForcing(1:end-1,1))~=0
    disp('timestamp of forcing data is not in regular intervals -> check, fix and restart')
    forcing.STATUS=0;
    return
else
    forcing.STATUS=1;
end

%here, consistency checks, RH->q calculation, set threhsolds for wind, etc. could be placed

% forcing.DATA.wind(forcing.DATA.wind<0.5)=0.5; %set min wind speed to 0.5 m/sec to avoid breakdown of turbulence
forcing.DATA.Lin(find(forcing.DATA.Lin==0)) = 5.67e-8 .* (forcing.DATA.Tair(find(forcing.DATA.Lin==0))+273.15).^4;

%set pressure to mean pressure at corresponding altitude (international
%altitude formula) if now provided by the forcing time series
if ~isfield(temp.FORCING.data, 'p')
    forcing.DATA.p=forcing.DATA.Tair.*0 + 1013.25.*100.*(1-0.0065./288.15.*forcing.PARA.altitude).^5.255;
else
    forcing.DATA.p = temp.FORCING.data.p;
end

if isempty(forcing.PARA.start_time) || ~ischar(forcing.PARA.start_time)
    forcing.PARA.start_time = forcing.DATA.timeForcing(1,1);
else
    forcing.PARA.start_time = datenum(forcing.PARA.start_time, 'dd.mm.yyyy');
end
if isempty(forcing.PARA.end_time) || ~ischar(forcing.PARA.end_time)
    forcing.PARA.end_time = floor(forcing.DATA.timeForcing(end,1));
else
    forcing.PARA.end_time = datenum(forcing.PARA.end_time, 'dd.mm.yyyy');
end


%initialize TEMP
forcing.TEMP.snowfall=0;
forcing.TEMP.rainfall=0;
forcing.TEMP.Lin=0;
forcing.TEMP.Sin=0;
forcing.TEMP.Tair=0;
forcing.TEMP.wind=0;
forcing.TEMP.q=0;
forcing.TEMP.p=0;
forcing.TEMP.Sin_dif = 0;
forcing.TEMP.Sin_dir = 0;
forcing.TEMP.snow_reservoir = 0;
