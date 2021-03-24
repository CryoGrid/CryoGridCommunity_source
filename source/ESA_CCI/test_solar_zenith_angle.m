Latitude = 70; %[70:-0.01:60];
Longitude = 10;%[100:0.01:110];
date = datenum(2000,6,12,12,0,0);
TZone = 0;

%DayOfYear = pvl_date2doy(Year, Month, Day);

Longitude = -Longitude;

Year = str2num(datestr(date, 'yyyy'));

DayOfYear = floor(date - datenum(Year,1,1) + 1);

%DecHours = Hour + Minute./60 + Second./3600;

DecHours = (date - floor(date)) .* 24;


RadtoDeg = 180 / pi;
DegtoRad = pi / 180;

Abber = 20/3600;
LatR = Latitude * DegtoRad;
UnivDate = DayOfYear + floor((DecHours + TZone)/24);
UnivHr = mod((DecHours + TZone), 24);
Yr = Year-1900;
YrBegin = 365 * Yr + floor((Yr-1)/4)-0.5;
Ezero = YrBegin + UnivDate;
T = Ezero / 36525;
GMST0 = 6/24 +38/1440 + (45.836 + 8640184.542 * T + 0.0929 * T.^2)/86400;
GMST0 = 360 * (GMST0 - floor(GMST0));
GMSTi = mod(GMST0 + 360*(1.0027379093 * UnivHr / 24),360);
LocAST = mod((360 + GMSTi - Longitude), 360);


EpochDate = Ezero + UnivHr / 24;
T1 = EpochDate / 36525;
ObliquityR = DegtoRad * (23.452294 - 0.0130125 * T1 - 0.00000164 * T1.^2 ...
    + 0.000000503 * T1.^3);
MlPerigee = 281.22083 + 0.0000470684 * EpochDate + 0.000453 * T1 .^ 2 + ...
    0.000003 * T1 .^ 3;
MeanAnom = mod((358.47583 + 0.985600267 * EpochDate - 0.00015 * T1 .^ 2 - ... 
    0.000003 * T1 .^ 3), 360);
Eccen = 0.01675104 - 0.0000418 * T1 - 0.000000126 * T1 .^ 2;
EccenAnom = MeanAnom;
E=0;

while max(abs(EccenAnom - E)) > 0.0001
    E = EccenAnom;
    EccenAnom = MeanAnom + RadtoDeg .* Eccen .* sin(DegtoRad .* E);
end

TrueAnom = 2 * mod(RadtoDeg * atan2(((1 + Eccen) ./ (1 - Eccen)).^ 0.5 .* tan(DegtoRad * EccenAnom / 2), 1), 360) ;
EcLon = mod(MlPerigee + TrueAnom, 360) - Abber ;
EcLonR = DegtoRad * EcLon;
DecR = asin(sin(ObliquityR) .* sin(EcLonR));
Dec = RadtoDeg * DecR;
RtAscen = RadtoDeg * atan2(cos(ObliquityR).*(sin(EcLonR)),cos(EcLonR));
HrAngle = LocAST - RtAscen ;
HrAngleR = DegtoRad .* HrAngle ; 

HrAngle = HrAngle - (360 .* sign(HrAngle) .* (abs(HrAngle) > 180));


SunAz = RadtoDeg .* atan2(-1 * sin(HrAngleR), cos(LatR) .* tan(DecR) - sin(LatR) .* cos(HrAngleR));
SunAz = SunAz + (SunAz < 0) * 360; %shift from range of [-180,180] to [0,360]
SunEl = asind(cos(LatR) .* cos(DecR) .* cos(HrAngleR) + sin(LatR) .* sin(DecR));

% Convert solar azimuth angle from [N,E,S,W]=[0,90,180,270] to [180, 90, 0
% -90], i.e. the same as the aspect and horizon angle system.
saz=deg2rad(SunAz); 
saz=(5*pi/2)-saz;
saz=saz-2.*pi.*(saz>2.*pi);
saz=saz+pi/2;
saz=saz-2.*pi.*(saz>pi);

% Calculate solar zenith angle from solar elevation angle
SunEl=deg2rad(SunEl);
szen=(pi/2)-SunEl;