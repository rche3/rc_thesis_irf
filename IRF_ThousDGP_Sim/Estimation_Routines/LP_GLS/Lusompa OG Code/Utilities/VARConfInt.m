
function [LO,HI,ME,algo] = VARConfInt(VAR,IRF_opt,ndraws,cnf_lvl,method,algo)

% Calculates confidence bands for two different significance levels
% for a VAR using bootstrapping.
%
% Usage:   [LO,HI,ME,algo] = VARConfInt(VAR,IRF_opt,...
%                                  ndraws,cnf_lvl,method,algo)
%
% Inputs:  VAR:     VAR results from VARreg
%          IRF_opt: Options of the IRFs (is the output of VAR)
% (opt.)   ndraws:  number of bootstrap draws (default: 100)
% (opt.)   cnf_lvl: Confidence level (default: 90)
% (opt.)   method:  'BS' for bootstrap (drawing from realised residuals) or
%                   'MC' for Monte Carlo (drawing from multivariate normal
%                   with var.-cov. matrix VAR.sigma). Default: 'BS'
% (opt.)   algo:    'F' for Fabio's, 'A' for Cesa-Ambrogio's.
%                   Note: 'F' algo is still on experimental stage; gives
%                         bad results!
% Outputs:
% Three arrays of dimension(h,v,s): horizon h, variable v, shock s
%          LO:     Lower confidence band
%          HI:     Upper confidence band
%          ME:     Median response
%
% by B. Kolb, Apr. 2015, based on code by A. Cesa Bianchi (VARirband.m)

% CAREFUL: I've changed the column order w.r.t. Ambrogio's codes!

%% Some renaming, checking and pre-allocating
hor      = IRF_opt.hor;
Ft       = VAR.Ft;  % rows are coefficients, columns are equations
nvars    = VAR.nvar;
nexog    = VAR.nexog;
nlags    = VAR.nlags;
determ   = VAR.determ;
nobs     = VAR.nobs;

if nexog~=0
    exog = VAR.X_EX;
end

% Defaults
if ~exist('ndraws','var')
    ndraws = 100;
end

if exist('cnf_lvl','var')
    pctl = cnf_lvl;
else
    pctl = 90;
end

if ~exist('method','var')
    method = 'BS';
end

if ~exist('algo','var')
    algo = 'A';
end

LO = NaN(hor,nvars,nvars);
HI = NaN(hor,nvars,nvars);
ME = NaN(hor,nvars,nvars);
data_draws = NaN(nobs,nvars); % forecasts of VAR plus innov terms
IRF = NaN(hor,nvars,nvars,ndraws); % IRF across draws

%% A) Standard method:
% Draw innovations, generate data, estimate VARs on them and get IRF and
% CIs as percentiles of draws
if strcmp(algo,'A') == 1
    wbar    = waitbar(0,'1','Name','Taking draws for confidence intervals');
    
    for tt=1:ndraws
        waitbar(tt/ndraws,wbar,sprintf('%12.9f',tt/ndraws))
        
        % draw random innovations
        if strcmp(method,'BS') == 1
            % BS: draw from realised residuals (with replacement):
            innov = VAR.resid(ceil(nobs*rand(nobs,1)),:);
        elseif strcmp(method,'MC') == 1
            % MC: draw from normal distribution
            innov = randn(nobs,nvars) * sqrtm(VAR.sigma);
        end
        
        % clear regressor vectors
        regr_lag = NaN(1,nvars*nlags); % lags for FC regressors
        if nexog~=0 % all FC regressors (incl. determ & exog.)
            regr_all = NaN(1,nvars*nlags+determ+nexog);
        else % all FC regressors (incl. determ)
            regr_all = NaN(1,nvars*nlags+determ);
        end
        
        % create artificial data set with random innovations
        for jj = 1:nobs
            
            if jj < nlags+1
                % a) Intialise first draw as data plus random innov
                data_draws(jj,:) = VAR.y(jj,:) + innov(jj,:);
                switch jj % create regressors recursively
                    case 1 % first iteration
                        regr_lag(1:nvars) = data_draws(jj,:);
                    otherwise % all others
                        regr_lag(1:jj*nvars) = ...
                            [data_draws(jj,:) ...
                            regr_lag(1:(jj-1)*nvars)];
                end
            else
                % b) Then, forecast from data and add innov
                for mm = 1:nvars
                    data_draws(jj,mm) = regr_all*Ft(1:end,mm) + innov(jj,mm);
                end
                % update regr_lag matrix
                regr_lag = [data_draws(jj,:) regr_lag(1,1:(nlags-1)*nvars)];
            end
            
            % Add determ. terms & exog. vars as regressors for forecast
            switch determ
                case 0
                    regr_all = regr_lag;
                case 1
                    regr_all = [1 regr_lag];
                case 2
                    T = (1:nobs)';
                    regr_all = [1 T(jj) regr_lag];
                case 3
                    T = (1:nobs)';
                    regr_all = [1 T(jj) T(jj).^2 regr_lag];
            end
            if nexog~=0
                regr_all(1,nvars*nlags+determ+1:end) = exog(jj,:);
            end
        end
        
        % get VAR draw estimated on created data (one per iteration tt)
        if nexog~=0
            VAR_draw = VARreg(data_draws(1:end,:),nlags,determ,exog);
        else
            VAR_draw = VARreg(data_draws(1:end,:),nlags,determ);
        end
        
        % calculate and save impulse responses for draw
        [IRF(:,:,:,tt), ~] = ...
            VARImpResp(VAR_draw,hor,IRF_opt.ident,IRF_opt.impact);
        
    end
    close(wbar);
end

%% B) Fabio's method
% Algorithm 4.4 in Canova (2007)
if strcmp(algo,'F') == 1
    np  = (nlags*nvars+determ)*nvars; % # parameters
    dof = nobs - np; % degrees of freedom
    
    B   = (VAR.x'*VAR.x)\VAR.x'*VAR.y;
    S   = VAR.sigma; % already corrected for d.o.f.!
    if strcmp(IRF_opt.ident,'short_run')
        A   = chol(S);
    else
        disp('This algorithm of IRF generation only works for short-run')
        disp('restrictions right now!');
        return
    end
    
    Sigma_draws = NaN(nvars,nvars,ndraws);
    Beta_draws = NaN((nvars*nlags+determ)*nvars,ndraws);
    
    for i=1:ndraws
        %BK: does dof really have to enter? Should be corrected for already
        PS = real(sqrtm(inv(S))); % principal square root of sigma
        u  = randn(nvars,dof); % random variation
        Sigma_draws(:,:,i) = inv(PS*(u*u')*PS');
        [uu, s, vv]        = svd((kron(Sigma_draws(:,:,i), inv(VAR.x'*VAR.x))));
        Beta_draws(:,i)    = reshape(B,np,1) + ...
            real(uu*sqrt(s)*vv')*randn(np,1);
    end
    % compute impulse responses
    basicimp = NaN(hor,nvars,nvars,ndraws);
    wimpu = NaN(nlags*nvars,nlags*nvars);
    
    % if strcmp(IRF_opt.ident,'short_run')
    %     C = A'\B';
    %     % get rid of determ. terms
    %     C = C(:,determ+1:end);
    % else
    %     disp('This method of IRF generation only works for short-run restrictions right now!');
    % end
    
    for d=1:ndraws
        %F: only draw beta; sigma is kept fixed to cut down on computations.
        bet = reshape(Beta_draws(:,d), nvars*nlags+determ, nvars);
        BB  = companionbk(bet,nlags,nvars,determ,'det_f','lag_blk');
        
        for j=1:hor
            % compute impulse responses
            if j==1
                wimpu(:,:,j) = eye(size(BB));
            else
                wimpu(:,:,j) = BB^(j-1);
            end
            basicimp(j,:,:,d) = A\wimpu(1:nvars,1:nvars,j); % * C
        end
    end
    
    IRF = sort(basicimp,4);
end


%% Compute percentiles (cnf_lvl upper and lower percentile + median)
pctl_inf = (100-pctl)/2;
pctl_sup = 100 - (100-pctl)/2;
LO(:,:,:) = prctile(IRF(:,:,:,:),pctl_inf,4);
HI(:,:,:) = prctile(IRF(:,:,:,:),pctl_sup,4);
ME(:,:,:) = prctile(IRF(:,:,:,:),50,4);

