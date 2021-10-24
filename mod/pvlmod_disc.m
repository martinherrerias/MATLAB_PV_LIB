function Kn = pvlmod_disc(Kt,AM)
% KN = PVLMOD_DISC(KT,AM) - Minimal implementation of PVL_DISC, to estimate beam-clearness-index
%   (KN = BNI/ENI) from clearness index (KT = GHI/EHI) and absolute air mass AM, using the DISC
%   model.
%
% REF:
% [1] Maxwell, E. L., "A Quasi-Physical Model for Converting Hourly Global Horizontal to Direct 
%   Normal Insolation", Technical Report No. SERI/TR-215-3087, Golden, CO: Solar Energy Research 
%   Institute, 1987.
% [2] J.W. "Fourier series representation of the position of the sun". Found at:
%   http://www.mail-archive.com/sundial@uni-koeln.de/msg01050.html on January 12, 2012
%
% See also: PVL_DISC, PVLMOD_DIRINT PVL_DATE2DOY PVL_EPHEMERIS PVL_ALT2PRES PVL_DIRINT PVL_LOUCHE 
%   PVL_ORGILL_HOLLANDS PVL_REINDL_1 PVL_REINDL_2 PVL_ERBS

    narginchk(2,2);
    compatiblesize(Kt,AM);
    
    % In the original function, Z is clipped to 87°, limiting AM < 15.2
    %   Ztemp(Z>87)=87; % Make sure you don't take the root of a negative.
    %   AM = 1./(cosd(Ztemp) + 0.15.*((93.885-Ztemp).^(-1.253))) .* pressure./ 101325;
    % Instead, clip AM to < argmin of Knc(AM) (below)
    AM0 = 17.2545;      
    AM(AM > AM0) = AM0;
    
    A = NaN(size(Kt));
    B = NaN(size(Kt));
    C = NaN(size(Kt));

    % For Kt > 0.6, set the components of A, B, and C
    i = Kt > 0.6;
    A(i) = -5.743 + Kt(i).*( 21.77 + Kt(i).*(-27.49 + Kt(i).*11.56));
   %B(i) =  41.40 + Kt(i).*(-118.5 + Kt(i).*( 66.05 + Kt(i).*31.90));
    B(i) =  41.40 + Kt(i).*(-118.5 + Kt(i).*( 66.06 + Kt(i).*31.90));
    C(i) = -47.01 + Kt(i).*( 184.2 + Kt(i).*(-222.0 + Kt(i).*73.81));

    % For Kt <= 0.6, set the components of A, B, and C
    i = Kt <= 0.6;
    A(i) = 0.512  + Kt(i).*( -1.56 + Kt(i).*( 2.286 - Kt(i)*2.222));
    B(i) = 0.370  + Kt(i).*( 0.962 );
    C(i) = -0.28  + Kt(i).*( 0.932 + Kt(i).*(-2.048 ));

    delKn = A + B.*exp(C.*AM);

    % The first term appears as "0.886" in some versions (see PVL_DISC) 
    Knc =  0.866 + AM.*(-0.122 + AM.*(0.0121 + AM.*(-0.000653 + AM.*0.000014)));
    Kn = Knc - delKn;
    Kn = max(0,Kn);
end