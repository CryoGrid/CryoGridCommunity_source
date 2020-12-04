function forcing = load_forcing_from_mat(forcing)

temp=load(['forcing/' forcing.PARA.filename], 'FORCING');

forcing.DATA.rainfall=temp.FORCING.data.rainfall.*forcing.PARA.rain_fraction;
forcing.DATA.snowfall=temp.FORCING.data.snowfall.*forcing.PARA.snow_fraction;
forcing.DATA.Tair = temp.FORCING.data.Tair;
forcing.DATA.Lin = temp.FORCING.data.Lin;
forcing.DATA.Sin = temp.FORCING.data.Sin;
forcing.DATA.q = temp.FORCING.data.q;
forcing.DATA.wind = temp.FORCING.data.wind;
forcing.DATA.timeForcing = temp.FORCING.data.t_span;


if std(forcing.DATA.timeForcing(2:end,1)-forcing.DATA.timeForcing(1:end-1,1))~=0
    disp('timestamp of forcing data is not in regular intervals -> check, fix and restart')
    forcing.STATUS=0;
    return
else
    forcing.STATUS=1;
end

%here, consistency checks, RH->q calculation, set threhsolds for wind, etc. could be placed

forcing.DATA.wind(forcing.DATA.wind<0.5)=0.5; %set min wind speed to 0.5 m/sec to avoid breakdown of turbulence
forcing.DATA.Lin(find(forcing.DATA.Lin==0)) = 5.67e-8 .* (forcing.DATA.Tair(find(forcing.DATA.Lin==0))+273.15).^4;

%set pressure to mean pressure at corresponding altitude (international
%altitude formula) if not provided by the forcing time series
if ~isfield(temp.FORCING.data, 'p')
    forcing.DATA.p=forcing.DATA.Tair.*0 + 1013.25.*100.*(1-0.0065./288.15.*forcing.PARA.altitude).^5.255;
else
    forcing.DATA.p = temp.FORCING.data.p;
end

if ~isfield(temp.FORCING.data, 'saltConcForcing')  %allows seeting the upper boundary of salt_c as a parameter
    forcing.DATA.saltConcForcing=forcing.DATA.Tair.*0 + forcing.PARA.saltConcForcing;
else
    forcing.DATA.saltConcForcing = temp.FORCING.data.saltConcForcing;
end


%initialize TEMP
forcing.TEMP.snowfall=0;
forcing.TEMP.rainfall=0;
forcing.TEMP.Lin=0;
forcing.TEMP.Sin=0;
forcing.TEMP.Tair=0;
forcing.TEMP.wind=0;
forcing.TEMP.RH=0;
forcing.TEMP.q=0;
forcing.TEMP.p=0;
forcing.TEMP.saltConcForcing=0;
