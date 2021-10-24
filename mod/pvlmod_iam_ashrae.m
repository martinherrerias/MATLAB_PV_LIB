function IAM = pvlmod_iam_ashrae(bzero, theta)
% PVLMOD_IAM_ASHRAE(b,theta) - Determine the incidence angle modifier using the ASHRAE
%   transmission model.
%
% MODIFIED from PVLIB_IAM_ASHRAE
%   - no warnings for 0 < theta or theta > acosd(b/(1+b))
%   - simpler parsing
%
% Input Parameters:
%   b - A parameter to adjust the modifier as a function of angle of
%     incidence. Typical values are on the order of 0.05 [3].
%   theta - The angle of incidence between the module normal vector and the
%     sun-beam vector in degrees. Theta must be a numeric scalar or vector.
%     For any values of theta where abs(theta)>90, IAM is set to 0. For any
%     values of theta where -90 < theta < 0, theta is set to abs(theta) and
%     evaluated. A warning will be generated if any(theta<0 or theta>90).
%     For values of theta near 90 degrees, the ASHRAE model may be above 1
%     or less than 0 due to the discontinuity of secant(theta). IAM values
%     outside of [0,1] are set to 0 and a warning is generated.
% 
% Output Parameters:
%   IAM - The incident angle modifier calculated as 1-b*(sec(theta)-1) as
%     described in [2,3]. IAM is a column vector with the same number of 
%     elements as the largest input vector.
%
% References: see PVL_IAM_ASHRAE
% 
% See also 
%       PVL_IAM_ASHRAE CHECKIAM	PVL_GETAOI PVL_EPHEMERIS PVL_SPA PVL_IAM_PHYSICAL

    narginchk(2,2);
    validateattributes(bzero,{'numeric'},{'real','positive','nonzero'},'','bzero');
    validateattributes(theta,{'numeric'},{'real'},'','theta');

    IAM = 1-bzero.*(secd(theta) - 1);
    thmax = acosd(bzero/(bzero+1));
    IAM(abs(theta)>=thmax)=0;
end
