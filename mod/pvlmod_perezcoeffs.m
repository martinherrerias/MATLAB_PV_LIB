 function [F1,F2,s11,s22,rho,Se] = pvlmod_perezcoeffs(BNI,DHI,GNE,sunel,model,method,varargin)
% [F1,F2] = PVLMOD_PEREZCOEFFS(BNI,DHI,GNE,SUNEL)
% [F1,F2] = PVLMOD_PEREZCOEFFS(BNI,DHI,GNE,SUNEL,MODEL)
% [F1,F2,S11,S22,RHO,SE] = PVLMOD_PEREZCOEFFS(BNI,DHI,GNE,SUNEL,MDL,METHOD) (§)
% [F1,F2,S11,S22,RHO,SE] = PVLMOD_PEREZCOEFFS(BNI,DHI,GNE,SUNEL,FF,METHOD,COV,SE) (§)
%
%   Returns the anisotropy indices for Circumsolar and Horizon-Brightening irradiance components
%   according to one of the Perez model coefficient sets. The first two syntaxes are simple
%   adaptations of PVPMC's PV_LIB original pvl_perez(...) function, returning coefficients only.
%
%   Syntaxes marked with (§) dettach from the original model, to use smooth interpolation instead
%   of binning, and/or estimate model uncertainty.
%
%   See PVLMOD_PEREZ to calculate tilted irradiance using F1, F2
%
% INPUT:
%   BNI - array of Beam Normal Irradiances
%   DHI - array of Diffuse Horizontal Irradiances
%   GNE - array of Normal Extraterrestrial irradiances
%   SUNEL - array of apparent solar elevation angles (degrees)
%
%       NOTE: arrays can be any size or dimension, but singleton-expansion-compatible 
%
%   MODEL - character string, Perez coefficient set. Default is 'allsitescomposite1990', other
%       options are {'allsitescomposite1988','sandiacomposite1988','usacomposite1988','france1988',
%       'phoenix1988','elmonte1988','albuquerque1988','capecanaveral1988','albany1988'}
%
%       Alternatively, MODEL can be a 6x8 matrix of custom-fit coefficients (see PEREZ_FIT)
%
%   METHOD (optional, EXPERIMENTAL) - using 'linear', 'makima', etc. as oposed to 'bin' (default), 
%       defines an interpolation method to replace the hard binning on 'epsilon' of the original 
%       Perez et al. model. This reduces artifacts (discontinuities) in the solutions, although
%       it detaches from the canonical model implementation.
%
%   COV,SE - use 48x48 covariance matrix COV for the model coefficients, and 8-vector of standard
%       errors SE, as returned by PEREZ_FIT, to estimate model uncertainty.
%
% OUTPUT:
%
%   F1 - Circumsolar Anisotropy index, F1 = 0 for overcast, F1 = 1 for black sky
%   F2 - Horizon-Brightening Anisotropy index, F2 = 0, no brightening
%
%   S11,S22,RHO,SE - estimated uncertainties, uncertainty correlation, and standard error of 
%       the resulting coefficients F1, F2, so that for each point k:
%
%       cov(F1(k),F2(k)) = s[1,r;r,1]s' for s = [S11(k),S22(k)]' and r = RHO(k)
%
%     If COV and SE are not provided as inputs, these values are estimated from the ensemble 
%     of all but one(*) coefficient sets, without correction for over- or under-dispersion.
%
% See also: PVL_PEREZ, PVLMOD_PEREZ, PVLMOD_PEREZGEOM, PVLMOD_HAYDAVIES
%
% (*) Leaving 'osage1988' out! From Perez et al. (1988):
% "With the exception of Osage, data are representativeof most solar geometries [...]
% the very poor results obtained withthe Osage-based model in Albuquerque and Phoenix are simply
% due to the fact that Osage, unlike the two southwestern sites, included only very few high-
% epsilon events andmean, because of the least square fitting method used to derive the coefficients,
% the resolution achieved for such events is totally unsatisfactory and allows for important 
% distortions. 
% Likewise, the small performance deterioration caused in all SNLA sites by the Albany-based model
% may be explained, in part, by the higher latitude of this site and the corresponding lack of very
% low solar zenith angle events."

    MODELS = {'allsitescomposite1990','allsitescomposite1988','sandiacomposite1988',...
             'usacomposite1988','france1988','phoenix1988','elmonte1988','albuquerque1988',...
             'capecanaveral1988','albany1988'}; % 'osage1988'
         
    % Parse input, and complete simulation options with defaults
    narginchk(4,8);
    assert(all(cellfun(@isnumeric,{BNI,DHI,GNE,sunel})),'Not numeric args. BNI, DHI, GNE, sunel');
    [BNI,DHI,GNE,sunel] = compatiblesize(BNI,DHI,GNE,sunel);
    
    if nargin < 5 || isempty(model), model = 'allsitescomposite1990'; 
    elseif isnumeric(model)
    % Use custom (fitted) coefficient matrix F = [F1c,F2c]'
        validateattributes(model,{'numeric'},{'real','finite','2d','size',[8,6]},'','F');
        F1custom = model(:,1:3);
        F2custom = model(:,4:6);
        model = 'custom';
        MODELS{end+1} = 'custom';
    else
        model = parselist(model,MODELS);
    end
    if nargin < 6, method = 'bin'; end

    if numel(varargin) > 0 && ~isempty(varargin{1})
        validateattributes(varargin{1},{'numeric'},{'real','finite','2d','size',[48,48]},'','V');
        isposdef = @(x) issymmetric(x) && all(eig(x) > -eps(1));
        assert(isposdef(varargin{1}),'Variance matrix should be symetric positive semi-definite');
        V = varargin{1};
    else
        V = [];
    end
    if numel(varargin) > 1 && ~isempty(varargin{2})
        Se = varargin{2};
        validateattributes(Se,{'numeric'},{'real','finite','vector','nonnegative'},'','SE');
        assert(isscalar(Se) || numel(Se) == 8,'Expecting scalar or 8-vector SE');
        if isscalar(Se), Se = repmat(Se,1,8); end
    else
    % TODO: fit with more data! these are based on EUO_Eugene (2019)
        Se = [20,40,50,50,50,40,25,15]; 
    end
        
    AMr = zeros(size(sunel));
    AMr(:) = pvl_relativeairmass(max(0,90-sunel(:)));

    kappa = 1.041; 
    zenith = (90-sunel)*pi/180; 
    epsilon = 1 + (BNI./DHI)./(1+kappa.*zenith.^3);
    epsilon(DHI == 0) = 1;

    delta = DHI.*AMr./GNE;
    delta(DHI == 0 | sunel == 0) = 0;
   
    X0 = [1, 1.065, 1.23, 1.5, 1.95, 2.8, 4.5, 6.2, Inf]; % bin edges (epsilon)

    [epsilon,delta,zenith] = compatiblesize(epsilon,delta,zenith);
    sz = size(epsilon);
    if numel(sz) == 2 && sz(1) == 1 && sz(2) > 1
    % make sure x(filter) is [n,1]
       epsilon = epsilon';
       delta = delta';
       zenith = zenith'; 
    end
    
    filter = epsilon >= 1 & delta >= 0 & zenith <= pi/2;
    n = nnz(filter);
    if n == 0
        F1 = NaN(sz); F2 = NaN(sz);
        if nargout > 2
            s11 = NaN(sz); s22 = NaN(sz); rho = NaN(sz); Se = NaN;
        end
        return;
    end
    epsilon = epsilon(filter);
    a = [ones(n,1),delta(filter),zenith(filter)];
    clear delta zenith

    if ~any(strcmpi(MODELS,model)) || (nargout > 2 && isempty(V))
        F1s = zeros(n,8,numel(MODELS));
        F2s = zeros(n,8,numel(MODELS));
        for j = numel(MODELS):-1:1

            model = MODELS{j};
            [F1c,F2c] = PerezCoefficients(model); % [8x3]
            F1s(:,:,j) = max(0,a*F1c');   % F11 + F12.*delta + F13.*zenith [nx8]
            F2s(:,:,j) = a*F2c';          % F21 + F22.*delta + F23.*zenith
            
            FF(:,:,j) = [F1c,F2c]; 
        end
    end
    
    switch lower(model)
    case 'custom'
        F1b = max(0,a*F1custom');
        F2b = a*F2custom';
    case 'min'
        F1b = min(F1s,[],3);
        F2b = min(F2s,[],3);
    case 'max'
        F1b = max(F1s,[],3);
        F2b = max(F2s,[],3);
    case 'median'
        F1b = median(F1s,3);
        F2b = median(F2s,3);
    case 'mean'
        F1b = mean(F1s,3);
        F2b = mean(F2s,3);
    otherwise
        [F1c,F2c] = PerezCoefficients(model);
        F1b = max(0,a*F1c');
        F2b = a*F2c';
    end
    % F1,F2 are at this point [nx8] arrays of coefficients (for the 8 bins of EPSILON)
    
    if nargout > 2 
        if isempty(V)
            FF = permute(FF,[3,2,1]);  % [6xnx8]
            V = cell(1,8);
            for j = 1:8
                V{j} = cov(FF(:,:,j));
            end
            V = blkdiag(V{:});
        end
        % V is covariance for F = [F11,F12,F13,F21,F22,F23] | bin 1, [F11 ... F23] | b2 , ..
        % Group all terms for F1jF1k into V11, F2jF2k into V22, and F1jF2k into V12
        V11 = V((0:2)' + (1:6:48),(0:2)' + (1:6:48));
        V22 = V((3:5)' + (1:6:48),(3:5)' + (1:6:48));
        V12 = V((0:2)' + (1:6:48),(3:5)' + (1:6:48)); % each [24,24]
    end

    switch lower(method)
    case 'bin'
        ebin = discretize(epsilon,X0);
        idx = sub2ind([n,8],(1:n)',ebin);
        F1 = F1b(idx);
        F2 = F2b(idx);
        
        if nargout > 2      
            % Variance elements must be weighted by repmat([1,d,z]·[1,d,z]',8,8) 
            % weight' 
            [aTa,r,c] = xTx_triu_elements(a);
            r = (ebin - 1)*3 + r';
            c = (ebin - 1)*3 + c';
            idx = sub2ind([24,24],r,c);
            
            s11 = sqrt(dot(aTa,V11(idx),2));
            s22 = sqrt(dot(aTa,V22(idx),2));
            rho = dot(aTa,V12(idx),2)./(s11.*s22);

            if ~isscalar(Se), Se = Se(ebin); end
        end
    otherwise
        
        %   [f, gof] = fit( (1:8)',X0(1:8)'-1,'power1','startpoint',[1e-3,3]);
        %   Xc = round(feval(f,1.5:8.5)+1,3)
        Xc = [1.015;1.091;1.295;1.708;2.425;3.552;5.204;7.505]; % bin centers
            
        ebin = interpmatrix(epsilon,Xc,'extrap','nearest','method1d',method);
        ebin = ebin./sum(ebin,2);
        
        F1 = dot(ebin,F1b,2);
        F2 = dot(ebin,F2b,2);

        if nargout > 2
        % Estimate covariance of F1,F2 based on covariance V of F 
        % If F = [F11,F12,F13,F21,F22,F23]', then [F1,F2]' = [a',0; 0,a']*F, with a' = [1,del,z]
        % Each element of the covariance of [F1 F2]' is then a'a.*Qjk, where V = [Q11,Q12;Q21;Q22]

            
            [mult,r,c] = xTx_triu_elements(ones(1,24));
            idx = sub2ind([24,24],r,c)';
            
            ra = mod(r-1,3)+1;
            ca = mod(c-1,3)+1;
            rb = floor((r-1)/3)+1;
            cb = floor((c-1)/3)+1;
            
            [~,iu,ia] = unique([ra,ca],'rows');
            aTa = a(:,ra(iu)).*a(:,ca(iu));
            
            [~,iu,ib] = unique([rb,cb],'rows');
            bTb = ebin(:,rb(iu)).*ebin(:,cb(iu));
            
            W = aTa(:,ia).*bTb(:,ib).*mult;

            s11 = sqrt(W*V11(idx)');
            s22 = sqrt(W*V22(idx)');
            rho = W*V12(idx)'./(s11.*s22);
            
            if ~isscalar(Se), Se = sum(ebin.*Se(:)',2); end
        end
        
    end
    F1 = revertfilter(full(F1),filter);
    F2 = revertfilter(full(F2),filter);
    
    if nargout > 2
        s11 = revertfilter(s11,filter);
        s22 = revertfilter(s22,filter);
        rho = revertfilter(rho,filter);
        if ~isscalar(Se), Se = revertfilter(full(Se),filter); end
    end
end

function [xTx,r,c] = xTx_triu_elements(x)
% For a row vector x, return a list of elements of triu(x'x), with all off-diagonal elements
% multiplied by 2, so that for a symmetric matrix B, B·(x'x) = dot(B(r,c),xTx)
% For a matrix X, do the same for each row X(:,j).

    n = size(x,2);
    r = nonzeros(triu(repmat(1:n,n,1)'));
    c = nonzeros(triu(repmat(1:n,n,1)));
    xTx = x(:,r).*x(:,c).*(1 + (r ~= c)');
end