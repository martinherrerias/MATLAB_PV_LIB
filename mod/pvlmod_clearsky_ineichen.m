function [CSGHI,CSDNI,CSDHI] = pvlmod_clearsky_ineichen(Location,t,TL,sunel,AMa,ENI,corrected)
% [ClearSkyGHI,ClearSkyDNI,ClearSkyDHI] = PVLMOD_CLEARSKY_INEICHEN(LOC,T) - Estimate Clear-Sky
%   Irradiance values for location LOC, at timesteps T, using interpolated monthly values for 
%   Linke Turbidity i.e. TL = LINKETURBIDITY(LOC,T)
%
% [ClearSkyGHI,ClearSkyDNI,ClearSkyDHI] = PVLMOD_CLEARSKY_INEICHEN(LOC,T*,[TL,SUNEL,AMA,ENI]) 
%   provide explicit Linke-Turbidity Factors TL for LOC at each time T; and/or precalculated solar
%   elevation angle SUNEL (degrees), absolute air-mass AMA, and/or normal extraterrestrial 
%   irradiance ENI.
%   (*) If all TL, SUNEL, AMA, and ENI are provided, T can be omitted.
%
% Changes vs PVL_CLEARSKY_INEICHEN:
%   + Input parsing (flexibility in input, precalculated airmass and extra. irradiance...)
%   + Linke Turbidity interpolation taken to external LINKETURBIDITY
%   + Removed artifact for low elevation angles
%
% Input Parameters:
%   LOC - a struct with scalar fields 'latitude', 'longitude' (degrees N, E), and 'altitude' (MASL)
%   T - [Nt,1] vector of UTC DATENUM time-steps-
%   TL - (optional) [Nt,1] vector of (AM2) Linke-Turbidity values.
%   SUNEL - (optional) [Nt,1] apparent SOLAR-elevation angle (degrees)
%   AMA - (optional) [Nt,1] absolute air-mass
%   ENI - (optional) [Nt,1] normal extraterrestrial irradiance
%
% Output:   
%   ClearSkyGHI - the modeled global horizonal irradiance in W/m²
%   ClearSkyDNI - the modeled direct normal irradiance in W/m²
%   ClearSkyDHI - the calculated diffuse horizonal irradiance in W/m²
%
% See also: PVL_CLEARSKY_INEICHEN, LINKETURBIDITY, PVL_CLEARSKY_HAURWITZ, EFFECTIVESOLARPOSITION

    narginchk(2,7);
    parselocation(Location,'optional',{'TimeZone','name'});
    
    if nargin < 3, TL = []; end
    if nargin < 4, sunel = []; end
    if nargin < 5, AMa = []; end
    if nargin < 6, ENI = []; end
    if nargin < 7, corrected = true; end
    
    if ~isempty(t)
    % Get any/all missing time-dependent variables
        if isempty(sunel)
            [~,~,sunel] = pvlmod_ephemeris(t,Location);
        end
        if isempty(TL), TL = linketurbidity(Location,t); end   
        if isempty(AMa)
            Pres = pvl_alt2pres(Location.altitude);
            AMr = pvl_relativeairmass(90-sunel,'kastenyoung1989');
            AMa = AMr.*Pres/101325;
        end
        % if isempty(ENI), ENI = pvl_extraradiation(t-datenum(year(t),0,0)); end
        if isempty(ENI), ENI = solarposition.extraradiation(t,1367); end
    end
    compatiblesize(sunel,TL,AMa,ENI);
    
    cosZ = max(0,sind(sunel));

    fh1=exp(Location.altitude.*(-1/8000));
    fh2=exp(Location.altitude.*(-1/1250));

    cg1=(0.0000509.*Location.altitude+0.868);
    cg2=0.0000392.*Location.altitude+0.0387;
    
    if corrected
    % PVLMOD: see: https://github.com/pvlib/pvlib-python/issues/435
    %
    % The correction term 0.01.*(AMc).^(1.8) creates an artifact where ClearSkyGHI shoots up at 
    % low elevation angles (most notably for low turbidity and/or altitude).
    % To remove it without introducing bias, air-mass is clipped to the point where the 
    % exponential correction becomes problematic:
    
        AM_max = (cg2.*(fh1+fh2.*(TL-1))/0.018).^1.25;
        AMc = min(AM_max,AMa);

        % AM_max is a minimum for the combined exponential terms, which shoot up afterwards. i.e.
        % AM_max = argmin(f(x)), for f = exp(-cg2*x*(fh1 + fh2*(tl - 1)) + 0.01*x^1.8))
        % It starts at AM = 2.6 (for turbidity = 1 and altitude = 0), equivalent to ~22.5° solar 
        % elevation, and goes up to O(10) for higher turbidities and altitudes. For most cases, it
        % corresponds to solar elevation angles below the 4° cutting threshold of the CIE 1994 quality
        % control guidelines that Ineichen & Perez seem to have used, i.e. they probably didn't see
        % beyond it.
        % Clipping at the minimum makes the adjustment smooth (C1), and leaves the results of the
        % model untouched for most operating conditions.
    else
       AMc = AMa; 
    end

    CSGHI = cg1.*ENI.*cosZ.*exp(-cg2.*AMc.*(fh1+fh2.*(TL-1))+0.01.*(AMc).^(1.8));
    CSGHI(cosZ<=0)=0;

    b = 0.664 + 0.163 ./ fh1;
    BncI = b .* ENI .* exp(-0.09 .* AMa .* (TL-1));
    BncI(cosZ<=0)=0;

    % Take the minimum of BncI and the equation given 
    CSDNI = min(BncI, CSGHI .* (1-(0.1 - 0.2 .* exp(-TL)) ./ (0.1 + 0.882 ./ fh1)) ./ cosZ);
    CSDHI = CSGHI - CSDNI .* cosZ;
    CSDHI(cosZ<=0)=0;
end
