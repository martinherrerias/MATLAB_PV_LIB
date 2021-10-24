function [GTI,ISO,CS,HB,ALB,BTI,Se] = pvlmod_perez(surftilt,surfaz,DHI,BNI,ENI,sunel,sunaz,albedo,varargin)
% [GTI,ISO,CS,HB,ALB,BTI] = PVLMOD_PEREZ(SURFTILT,SURFAZ,DHI,BNI,ENI,SUNEL,SUNAZ,ALBEDO,MODEL)
%
%   Estimate non-shaded global tilted irradiance based one of the Perez models.
%   Adaptation of PVPMC's PV_LIB original PVL_PEREZ(...) function, allowing for input arrays of
%   any size or dimension (singleton-expansion-compatible)* with missing values
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
% TODO: DOCUMENT UNCERTAINTY!
%
% CHANGES WITH RESPECT TO PVL_PEREZ: aside from input-output format, and the fact albedo and beam
%   components are included, the only change with respect to PV_LIB's PVL_PEREZ function is the use
%   of constant air mass at horizon (37.0) for all points with SUNEL < 0.
%
% References
%   [1] Loutzenhiser P.G. et. al., 2007. Empirical validation of models to compute
%   solar irradiance on inclined surfaces for building energy simulation, 
%   Solar Energy vol. 81. pp. 254-267.
%   [2] Perez, R., Seals, R., Ineichen, P., Stewart, R., Menicucci, D., 1987. A new
%   simplified version of the Perez diffuse irradiance model for tilted
%   surfaces. Solar Energy 39 (3), 221–232.
%   [3] Perez, R., Ineichen, P., Seals, R., Michalsky, J., Stewart, R., 1990.
%   Modeling daylight availability and irradiance components from direct
%   and global irradiance. Solar Energy 44 (5), 271–289.
%   [4] Perez, R. et. al 1988. The Development and Verification of the
%   Perez Diffuse Radiation Model,.SAND88-7030, Sandia National
%   Laboratories.
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
    
    [a,~,ISO,CS,HB,ALB] = pvlmod_perezgeom(surftilt,surfaz,sunel,sunaz);
    
    % Check input that is not already parsed by PVLMOD_PEREZCOEFFS
    assert(isnumeric(albedo),'Non-numeric args');
    compatiblesize(a,albedo,F1);

    GHI = max(0,DHI + BNI.*sind(sunel));
    
    A = DHI.*(CS - ISO);
    B = DHI.*HB;
    
    ISO = (1-F1).*DHI.*ISO;   %	Isotropic component: (1-F1)·DHI·(1 + k'u)/2	
    CS = F1.*DHI.*CS;         % Circumsolar component: F1·DHI·a/b
    HB = F2.*DHI.*HB;         % Hor. brightening component: F2·DHI·sqrt(1-(k'u)²) 
    BTI = BNI.*a;             % Beam component
    ALB = albedo.*GHI.*ALB;   %	Albedo component: rho·GHI·(1 - k'u)/2
    
    GTI = ISO + CS + HB + ALB + BTI;

    if nargout > 6
    % Does not include covariance (several sensors at the same time)!
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
