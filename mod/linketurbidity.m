function TL = linketurbidity(Loc,t)
% TL = LINKETURBIDITY(LOC) - lookup monthly Linke Turbidity Factors for location LOC.
% TL = LINKETURBIDITY(LOC,T) - interpolated Factors for LOC on the dates specified by T.
%
%   The .mat file 'LinkeTurbidities.mat' must be available on MATLAB's path, and contain a single
%   2160 x 4320 x 12 uint8 matrix variable called 'LinkeTurbidity'. Rows represent global latitudes
%   from 90 to -90 degrees; columns represent global longitudes from -180 to 180; and depth (third
%   dimension) represents months of the year from January (1) to December (12).
% 
%   CHANGES from PVL_CLEARSKY_INEICHEN:
% 
%       LINKETURBIDITY uses cubic interpolation, first to calculate monthly averages based on
%       location, and then for daily values in-between monthly average-days.
%
% See also: PVLMOD_DETECT_CLEAR_TIMES, PVLMOD_CLEARSKY_INEICHEN, LINKETURBIDITIES.mat, 
%   EFFECTIVESOLARPOSITION

    narginchk(1,2);
    parselocation(Loc,'optional',{'TimeZone','name'});
    if nargin < 2 || isempty(t)
    % Average days for months (*). 2019 is hip, but any non-leap year should do. 
        DoY = solarposition.mdoy();
    else
        DoY = solarposition.doy(t);
    end
    assert(~isempty(which('LinkeTurbidities.mat')),'LinkeTurbidities.mat not found.');

    % Create and evaluate cubic interpolant for monthly averages at the given location
    load('LinkeTurbidities.mat','LinkeTurbidity');
    LinkeTurbidity = double(flip(LinkeTurbidity,1))/20; % LT/20, switch 90:-90 to -90:90
    x = -90+1/24:1/12:90-1/24;
    y = -180+1/24:1/12:180-1/24;
    F = griddedInterpolant({x,y,1:12},LinkeTurbidity,'linear');
    TL = F(ones(1,12)*Loc.latitude,ones(1,12)*Loc.longitude,1:12);    % 12 vector

    % Interpolate in time from the 12 monthly averages, using average days for Months (*)
    % (as recommended from Klein, 1977 - Duffie Beckman 3rd Ed.)
    % Include nov-dec for the year before, and jan-feb for the year after, so that f(0) = f(365)
    
    n0 = solarposition.mdoy([11:12,1:12,1:2]) + 365.*[-1,-1,zeros(1,12),1,1];
    TL = TL([11,12,1:12,1,2]);

    TL = interp1(n0,TL,DoY,'pchip');
end
