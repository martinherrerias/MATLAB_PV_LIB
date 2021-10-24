function [SunAz, SunEl, ApparentSunEl, SolarTime, Dec]= pvlmod_ephemeris(Time, Location, varargin)
% [SunAz, SunEl, ApparentSunEl, SolarTime]=pvlmod_ephemeris(Time, Location)
% [SunAz, SunEl, ApparentSunEl, SolarTime]=pvlmod_ephemeris(Time, Location, Pressure)
% [SunAz, SunEl, ApparentSunEl, SolarTime]=pvlmod_ephemeris(Time, Location, Pressure, Temperature)
% [SunAz, SunEl, ApparentSunEl, SolarTime]=pvlmod_ephemeris(Time, Location, 'temperature', Temperature)
%
% Calculate the position of the sun given time, location, and optionally pressure and temperature
%
% Description: behaviour is the same as PVL_EPHEMERIS, except for the following changes
%
% PVLMOD: changes with respect to PVL_EPHEMERIS
%   + output declination
%   + use altitude-corrected mean pressure, if missing (with warning)
%   + warning if using 12° temperature as default
%   + looser-form, stricter-content parsing rules (allow NaNs)
%   + flexible Time format
%
% Input Parameters:
%   Time - see PARSETIME (PVLMOD)
% 
%   Location is a struct with the following elements 
%   Location.latitude = vector or scalar latitude in decimal degrees (positive is
%     northern hemisphere)
%   Location.longitude = vector or scalar longitude in decimal degrees (positive is 
%     east of prime meridian)
%   Location.altitude = an optional, used if pressure is not provided.
% 
% Output Parameters:
%   SunAz = Azimuth of the sun in decimal degrees from North. 0 = North to 270 = West
%   SunEl = Actual elevation (not accounting for refraction)of the sun 
%     in decimal degrees, 0 = on horizon. The complement of the True Zenith
%     Angle.
%   ApparentSunEl = Apparent sun elevation accounting for atmospheric 
%     refraction. This is the complement of the Apparent Zenith Angle.
%   SolarTime = solar time in decimal hours (solar noon is 12.00).
%
% References
%   Grover Hughes' class and related class materials on Engineering 
%   Astronomy at Sandia National Laboratories, 1985.
%
% See also PVL_MAKETIMESTRUCT PVL_MAKELOCATIONSTRUCT PVL_ALT2PRES
%          PVL_GETAOI PVL_SPA

narginchk(2,4);
Time = parsetime(Time);

p = inputParser;
p.addOptional('pressure',[], @(x) isnumeric(x) && isreal(x) && ~any(x < 53000 | x > 114000,'all'));
p.addOptional('temperature',12, @(x) isnumeric(x) && isreal(x) && ~any(x < -90 | x > 55,'all'));
p.parse(varargin{:});
pressure = p.Results.pressure; 
temperature = p.Results.temperature;

if nargout > 2
    if isempty(pressure)
        foo = naptime('parselocation:timezone'); %#ok<NASGU>
        Location = parselocation(Location,'-soft'); % assume altitude = 0 with warning, if missing
        pressure = pvl_alt2pres(Location.altitude);
    else
        Location = parselocation(Location);
    end
    if isempty(temperature)
        temperature = 12;
        warning('spa:temperature','Using 12°C for diffraction correction');
    end
    compatiblesize(pressure,temperature,Time);
end

TZone = hours(tzoffset(Time)); 
DayOfYear = floor(days(Time -dateshift(Time,'start','year'))) + 1;
% DayOfYear = pvl_date2doy(year(Time),month(Time),day(Time));
DecHours = hour(Time) + minute(Time)./60 + second(Time)./3600;

RadtoDeg = 180 / pi;
DegtoRad = pi / 180;

Abber = 20/3600;
LatR = Location.latitude * DegtoRad;
UnivDate = DayOfYear + floor((DecHours - TZone)/24);
UnivHr = mod((DecHours - TZone), 24);

Yr = year(Time)-1900;
YrBegin = 365 * Yr + floor((Yr-1)/4)-0.5;
Ezero = YrBegin + UnivDate;
T = Ezero / 36525;
GMST0 = 6/24 +38/1440 + (45.836 + 8640184.542 * T + 0.0929 * T.^2)/86400;
GMST0 = 360 * (GMST0 - floor(GMST0));
GMSTi = mod(GMST0 + 360*(1.0027379093 * UnivHr / 24),360);
LocAST = mod((360 + GMSTi + Location.longitude), 360);
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

SolarTime = (180 + HrAngle) / 15;

if nargout > 2
  
    % Calculate the refraction of the sun until the actual center of the sun is
    % 1 degree below the horizon.
    TanEl = tan(DegtoRad * SunEl);
    
    Refract = zeros(size(SunEl));
    f = SunEl > 5 & SunEl <= 85;
    Refract(f) = 58.1 ./ TanEl(f) - 0.07*TanEl(f).^(-3) + .000086*TanEl(f).^(-5);   
    f = SunEl > -0.575 & SunEl <=5;
    Refract(f) = SunEl(f) .* (-518.2 + SunEl(f) .* (103.4 + SunEl(f) .* (-12.79 + SunEl(f) .* 0.711))) +1735; 
    f = SunEl > -1 & SunEl <= -0.575;
    Refract(f) = -20.774 ./ TanEl(f);

    Refract = Refract .* (283 ./ (273 + temperature)) .* pressure ./ 101325 ./ 3600;

    % Generate apparent sun elevation including refraction
    ApparentSunEl = SunEl + Refract;
end



    

    