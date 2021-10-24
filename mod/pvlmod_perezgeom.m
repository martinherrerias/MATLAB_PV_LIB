function [a,b,ISO,CS,HB,ALB] = pvlmod_perezgeom(surftilt,surfaz,sunel,sunaz)
% [a,b,ISO,CS,HB,ALB] = PVLMOD_PEREZGEOM(SURFTILT,SURFAZ,SUNEL,SUNAZ)

% See also: PVLMOD_PEREZ, PVLMOD_PEREZCOEFFS

    % Parse input, and complete simulation options with defaults
    narginchk(4,4);
    
    % Check input that is not already parsed by PVLMOD_PEREZCOEFFS
    assert(all(cellfun(@isnumeric,{surftilt,surfaz,sunel,sunaz})),'Non-numeric args');
    compatiblesize(surftilt,surfaz,sunel,sunaz);
    
    % Allow dot-product-like behaviour for arrays of any (compatible) size, without creating
    % or processing unnecessary data copies:
    sph2cartC = @(az,el) {cosd(az).*cosd(el),sind(az).*cosd(el),sind(el)};
    dotC3 = @(A,B) A{1}.*B{1} + A{2}.*B{2} + A{3}.*B{3};
    
    % U contains unit vectors pointing away from the surface
    % S contains unit vectors pointing towards the sun
    S = sph2cartC(sunaz,sunel);
    U = sph2cartC(surfaz,90-surftilt);
    
    a = max(0,dotC3(S,U));
    b = max(0.087,S{3}); 
    
    ISO = (1 + U{3})/2; 
    CS = a./b;
    HB = hypot(U{1},U{2});
    ALB = 0.5*(1 - U{3});
    
    % Temps & Coulson
    % z = 90-sunel;
    % ALB = 0.5*albedo.*GHI.*(1 - U{3}).*(1+sind(z/2).^2).*abs(cosd(surfaz-sunaz));
end
