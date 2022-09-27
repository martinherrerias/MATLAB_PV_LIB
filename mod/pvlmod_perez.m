function [GTI,ISO,CS,HB,ALB,BTI,Se] = pvlmod_perez(surftilt,surfaz,DHI,BNI,ENI,sunel,sunaz,albedo,varargin)
% [GTI,ISO,CS,HB,ALB,BTI] = PVLMOD_PEREZ(SURFTILT,SURFAZ,DHI,BNI,ENI,SUNEL,SUNAZ,ALBEDO,MODEL)
%
%   Adaptation of PVPMC's PV_LIB original PVL_PEREZ(...) function, allowing for input arrays of
%   any size or dimension (singleton-expansion-compatible)* with missing values
%
% [GTI,ISO,CS,HB,ALB,BTI,SE] = PVLMOD_PEREZ(...) - estimate model uncertainty SE from the ensemble
%   of coefficient sets.
%
% [GTI,ISO,CS,HB,ALB,BTI,SE] = PVLMOD_PEREZ(...,MODEL,METHOD,COV,SE) - Allows the use of a custom
%   coefficient set (passed as MODEL), with known covariance and standard error. See PEREZ_FIT and
%   PVLMOD_PEREZCOEFFS.
%
%   SURFTILT - Surface tilt (elevation) angles (degrees)
%   SURFAZ - Surface azimuth angles (degrees) any convention, consistent with SUNAZ (§)
%	DHI - Diffuse Horizontal Irradiances
%	BNI - Beam Normal Irradiances
%	ENI - Extraterrestrial Normal irradiances
%	SUNEL - effective apparent solar elevation angles (degrees)
%	SUNAZ - effective apparent solar azimuth angles (degrees) convention consistent with SURFAZ (§)
%   ALBEDO - ground reflectivity (default 0.2).
%
%   MODEL - character string, Perez coefficient set. Default is 'allsitescomposite1990', other
%   options are {'allsitescomposite1988','sandiacomposite1988','usacomposite1988','france1988',
%   'phoenix1988','elmonte1988','osage1988','albuquerque1988','capecanaveral1988','albany1988'}
%   Alternatively, MODEL can be a 6x8 matrix of custom-fit coefficients (see PEREZ_FIT)
%
%   METHOD (optional, EXPERIMENTAL) - using 'linear', 'makima', etc. as oposed to 'bin' (default), 
%       defines an interpolation method to replace the hard binning on 'epsilon' of the original 
%       Perez et al. model. This reduces artifacts (discontinuities) in the solutions, although
%       it detaches from the canonical model implementation.
%
%   COV,SE - use 48x48 covariance matrix COV for the model coefficients, and 8-vector of standard
%       errors SE, as returned by PEREZ_FIT, to estimate model uncertainty.
%
%  GTI - Global (unshaded) tilted irradiance
%  ISO, CS, HB - Isotropic, Circumsolar, Horizon-Brightening components
%  ALB - Albedo component
%  BTI - Beam tilted irradiance
%
%   (*) NOTE: the size of the output arrays is the singleton-expanded combination of the inputs,
%       namely S = COMPATIBLESIZE(SURFTILT,...,ALBEDO,'-size').
%       e.g. if DHI,BNI,... are [N,1] vectors and SURFEL, SURFAZ are [1,M] vectors, results will
%       be [N,M] matrices. 
%
% See also: PVL_PEREZ, PVLMOD_PEREZCOEFFS, PVLMOD_HAYDAVIES, POAIRRADIANCE, EFFECTIVESOLARPOSITION

    % Parse input, and complete simulation options with defaults
    narginchk(7,12);
    if nargin < 8 || isempty(albedo)
        warning('Assuming 0.2 albedo');
        albedo = 0.2; 
    end
    
    if nargout > 6
        [F1,F2,sF1,sF2,rF1F2,Se] = pvlmod_perezcoeffs(BNI,DHI,ENI,sunel,varargin{:});
    else
        [F1,F2] = pvlmod_perezcoeffs(BNI,DHI,ENI,sunel,varargin{:});
    end
    
    [a,~,Fiso,Fcs,Fhb,Falb] = pvlmod_perezgeom(surftilt,surfaz,sunel,sunaz);
    
    % Check input that is not already parsed by PVLMOD_PEREZCOEFFS
    assert(isnumeric(albedo),'Non-numeric args');
    compatiblesize(a,albedo,F1);

    GHI = max(0,DHI + BNI.*sind(sunel));
    
    ISO = (1-F1).*DHI.*Fiso;   %	Isotropic component: (1-F1)·DHI·(1 + k'u)/2	
    CS = F1.*DHI.*Fcs;         % Circumsolar component: F1·DHI·a/b
    HB = F2.*DHI.*Fhb;         % Hor. brightening component: F2·DHI·sqrt(1-(k'u)²) 
    BTI = BNI.*a;              % Beam component
    ALB = albedo.*GHI.*Falb;   %	Albedo component: rho·GHI·(1 - k'u)/2
    
    GTI = ISO + CS + HB + ALB + BTI;

    if nargout > 6
    % Does not include covariance (several sensors at the same time)!
    
        A = DHI.*(Fcs - Fiso);
        B = (DHI.*Fhb);
    
        Se = Se.*sind(surftilt);
        Se = sqrt(Se.^2+(sF1.*A).^2+(sF2.*B).^2 + 2*rF1F2.*A.*B.*sF1.*sF2);
    end
end

% % TEST
% [MD,t,Loc] = getMeteoData('*.meteo');
% [MD,t,SP,Loc] = completemeteodata(MD,t,Loc);
% 
% f = ~(MD.flags.GHI > 0 | MD.flags.DHI > 0 | MD.flags.BNI > 0);
% surftilt = 60;
% surfaz = 135;
% SP.Az(SP.Az < 0) = SP.Az(SP.Az < 0) + 360;
% AMr = pvl_relativeairmass(90-SP.El(f));
% [~,ISO,CS,HB] = pvl_perez(surftilt,surfaz,MD.DHI(f),MD.BNI(f),MD.ENI(f),90-SP.El(f),SP.Az(f),AMr);
% [~,ISOm,CSm,HBm] = pvlmod_perez(surftilt,surfaz,MD.DHI(f),MD.BNI(f),MD.ENI(f),SP.El(f),SP.Az(f),0.2);
% plot([ISO,CS,HB],[ISOm,CSm,HBm],'.')

