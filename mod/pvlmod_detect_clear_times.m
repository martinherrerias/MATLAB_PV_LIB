function [CS,alpha] = pvlmod_detect_clear_times(GHI,CSGHI0,varargin)
% [CS,alpha] = PVLMOD_DETECT_CLEAR_TIMES(GHI,CSGHI,DT) 
% [CS,alpha] = PVLMOD_DETECT_CLEAR_TIMES(...,'name',Value)
%
%   Identify time samples with irradiance  GHI consistent with clear-sky CSGHI in a uniform time-
%   series with interval DT (MINUTES).
%
%   PVLMOD: changes vs PVL_DETECT_CLEAR_TIMES
%       + Clear sky model provided as input, e.g. from PVLMOD_CLEARSKY_INEICHEN
%       + Uniform time series only (performance)
%       + Name,Value - optional parameters for all thresholds
%       + EXPERIMENTAL: Scaling of slopes and window-length for dt ~= 1 min
%
% Inputs:   
%   GHI - Nx1 vector of GHI values (W/m²)
%   CSGHI - Nx1 vector of Clear-Sky GHI estimates (W/m²), e.g. from PVLMOD_CLEARSKY_INEICHEN
%   LOC - struct with the following elements, note that all
%      elements must be scalars in this application.
%   LOC.latitude - scalar latitude in decimal degrees (positive is
%      northern hemisphere)
%   LOC.longitude - scalar longitude in decimal degrees (positive is
%      east of prime meridian)
%   LOC.altitude - scalar height above sea level in meters.
%      While altitude is optional in many uses, it is required in this
%      model implementation.
%
% Output:   
%   CS - column vector of logical values (True = clear sample)
%   alpha - scaling factor applied to Ineichen model to minimize the error with CSGHI
%   
% PVLMOD_DETECT_CLEAR_TIMES(...,Name,Value) - Override default threshold values for each criterion, 
%   and/or set sliding-time-window length / nominal time-step:
%
%   'WinLength' - length of sliding time window (minutes), integer.
%       PVLMOD: Default is round(10·dt^(-0.3))·dt, where dt is the nominal time-step (see below).
%       10·t^(-0.3) yields ~10 time-steps for 1-min dt, ~3 time-steps for 60 min dt.
%                                        
%   'TimeStep' - nominal time step (minutes) between each GHI sample, integer. 
%
%   'MeanDiff': threshold value in W/m² for first criterion, i.e., agreement between  
%       mean values of GHI in each interval, see Eq. 6 in [1]
%       Default (recommended) is 75 W/m² for 1-min samples, 10-min window
%
%   'MaxDiff': threshold value in W/m² for second criterion, i.e., agreement between
%       maxima of GHI values in each interval, see Eq. 7 in [1]
%       Default (recommended) is 75 W/m² for 1-min samples, 10-min window
%
%   'LineLengthLower' and 'LineLengthUpper': threshold values (unitless) for third criterion
%       i.e. agreement between line lengths (LL) of GHI data and clear sky GHI, Eq. 8 in [1].
%       Criterion is satisfied when LineLengthLower <= LL <= LineLengthUpper
%       Default (recommended) is -5 <= LL <= 10 using LL² = (dP/[W])² + (dt/[min])²  See note [§]
%
%   'VarDiff': threshold value in 1/min [§] for the fourth criterion, i.e. agreement between
%       normalized standard deviations of rate of change in irradiance, Eqs. 9-11 in [1]
%       Default (recommended) is 0.005 1/min for 1-min time-steps and a 10-min window
%
%   'SlopeDev': threshold value in W/m² for the fifth criterion, i.e. agreement between largest
%       magnitude of change in successive GHI values, see Eq. 12 through 14 in [1]
%       Default (recommended) is 8 W/m² for 1-min samples, 10-min window
%
%   'TolAlpha': 1e-4, stop criteria for iterative search (alpha = scale factor GHI(CS)\CSGHI)
%   'MaxIter': 20, maximum number of iterations for alpha search.
%     
% References
%   [1] Reno, M.J. and C.W. Hansen, "Identification of periods of clear sky
%   irradiance in time series of GHI measurements" Renewable Energy, v90, 
%   p. 520-531, 2016.
%
% Notes:
%   Initial implementation by Matthew Reno. Modifications for computational efficiency by 
%   Joshua Patrick and Curtis Martin.
%
%   [§] PVLMOD: it's not clear from the paper what the units are, but in any case it seems 
%   reasonable to assume that they're 1/min (the standard time-step), not 1/seconds, as documented.
%
% TODO: Check Jamie Bright's https://github.com/JamieMBright/csd-library
%
% See also: PVL_DETECT_CLEAR_TIMES, PVLMOD_CLEARSKY_INEICHEN, FITCLEARSKY

    narginchk(3,Inf);

    OPT.TimeStep = [];
    OPT.MeanDiff = 75;
    OPT.MaxDiff = 75;
    OPT.LineLengthLower = -5; 
    OPT.LineLengthUpper = 10;
    OPT.VarDiff = 0.005;
    OPT.SlopeDev = 8;

    OPT.WinLength = [];    
    OPT.TolAlpha = 1e-4;
    OPT.MaxIter = 20;

    OPT = getpairedoptions(varargin,OPT,'dealrest',1);
    dt = OPT.TimeStep;
    if isduration(dt), dt = minutes(dt); end
    validateattributes(dt,{'numeric'},{'scalar','real','finite','positive','<=',60},'','TimeStep');
    
    % PVLMOD: adjust window length to timestep, 10 samples per window @ 1 min, 3 samples @ 60 min
    if isempty(OPT.WinLength), OPT.WinLength = round(10*dt^(-0.3))*dt; end    
    assert(mod(OPT.WinLength, dt) == 0,'OPT.WinLength must be an integer multiple of dt');
    if dt ~= 1
        warning('detectcleartimes:notoneminute',...
            'Using Reno & Hansen (2016) clear-sky detection for dt ~= 1 min'); 
    end
    
    Nt = numel(GHI);
    assert(numel(CSGHI0) == Nt,'Inconsistent GHI and CSGHI');

    % Create index matrix for windows
    samples_per_window = OPT.WinLength / dt;
    H = hankel(1:samples_per_window,samples_per_window:Nt)';
    clearMean = mean(CSGHI0(H),2);
    clearMax = max(CSGHI0(H),[],2);
    %clearSlope = diff(csGHI0(H),1,2);
    clearSlope = diff(CSGHI0(H),1,2)/dt; % PVLMOD: added dt-scaling [§]
    clearMaxSlope = max(abs(clearSlope),[],2);

    % Calculate parameters for measured GHI
    measuredMean = mean(GHI(H),2);
    measuredMax = max(GHI(H),[],2);
    %measuredSlope = diff(GHI(H),1,2);
    measuredSlope = diff(GHI(H),1,2)/dt; % PVLMOD: added dt-scaling [§]

    measuredMaxSlope = max(abs(measuredSlope),[],2);
    measuredLineLength = sum(sqrt(measuredSlope.^2 + 1)*dt,2); % [§]

    % Normalized variance of measured GHI < Vardeiff criterion doesn't depend on
    % clear sky GHI, compute outside of iteration on CS model adjustment alpha
    c4 = std(measuredSlope,[],2) ./ measuredMean < OPT.VarDiff;
    
    alpha = 1; % initialize
    c0 = measuredMean > 0 & clearMean > 0;

    % Start Iterative process of finding clear days, fitting clear sky model, finding clear days, ....
    % Continue iterations until not finding any more clear days or 20 iterations
    for ii = 1:OPT.MaxIter

        % Adjust Clear Sky Model parameters for scaling by alpha
        %clearLineLength = sum(sqrt((alpha*clearSlope).^2 + (dt*ones(size(measuredSlope))).^2 ));
        clearLineLength = sum(sqrt((alpha.*clearSlope).^2 + 1)*dt,2);  %[§]

        c1 = abs(measuredMean - alpha.*clearMean) < OPT.MeanDiff;
        c2 = abs(measuredMax - alpha.*clearMax) < OPT.MaxDiff;
        c3 = (OPT.LineLengthLower < measuredLineLength - clearLineLength) & (measuredLineLength - clearLineLength < OPT.LineLengthUpper) ;
        c5 = (measuredMaxSlope - alpha.*clearMaxSlope) < OPT.SlopeDev;
        c6 = clearMean ~= 0 & ~isnan(clearMean);            % window includes some daylight

        % Identify windows meeting all five (six) criteria
        clearWindows = c0 & c1 & c2 & c3 & c4 & c5 & c6;

        % Daily clearness is proportion of samples that are clear
        % dailyClearness = histcounts(Time(clearWindows),days) ./ histcounts(Time(csGHI>25),days);

        % Save current state
        previousAlpha = alpha;

        % Get new alpha by adjusting clear sky model to match clear times
        % Need to minimize RMSE between GHI and alpha * csGHI0
        CS = ismember((1:Nt)',H(clearWindows,:));
        
        RMSE = @(x) std(GHI(CS) - x*CSGHI0(CS),'omitnan');
        if isnan(RMSE(alpha)), break; end % e.g. ~any(CS)
        alpha = fminsearch(RMSE, alpha);  % update scaling factor on clear sky model

        if abs(alpha - previousAlpha) < OPT.TolAlpha, break; end

    end  % for loop
end  % function
