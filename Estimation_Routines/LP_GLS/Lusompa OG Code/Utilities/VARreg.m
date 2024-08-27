
function VARresults = VARreg(data,nlags,determ,varx,nlagx)

% Estimates a vector autogressive (VAR) system by OLS.
%
% Usage:   VARresults = runVAR(data,nlags,determ,varx,nlagx)
%
% Inputs:  data:   Matrix of data (assumed: # periods > # variables)
%          nlags:  Number of lags
%          determ: (opt.) Deterministic regressors: 0=none, 1=const,
%                  2=const&trend, 3=const&trend^2 (default is 1)
%	       varx:   (opt.) Matrix of exogenous variables (T x ncoex)
%          nlagx:  Number of lags for exogenous variables (default = 0)
%
% Output:  Structure "VARresults" containing the fields
% ~.nobs:       # periods (without lags)
% ~.nvar:       # endogenous variables
% ~.nexog:      # exogenous variables
% ~.nlags:      # lags for endogenous variables
% ~.nlags_exo:  # lags for exogenous variables
% ~.ncoeff:     # coeff's to estimate per eq. (nlags*nvar)
% ~.ncoeff_exo: # exogenous variables
% ~.ncoeff_tot: total # coeff's to estimate per equation
%                        (nlags*nvar + ncoeff_exo + determ)
% ~.determ:     0 none; 1 const; 2 const&trend; 3 const&trend^2
% ~.eqXX:       OLS results for each equation
% ~.Ft:         Matrix of estimated coeff's (as in VARform.m)
% ~.sigma:      VCV matrix of residuals (adjusted for dof)
% ~.resid:      Matrix of reduced form residuals
% ~.x:          Independent variable (as in VARform.m)
% ~.y:          Dependent variable
% ~.x_exo:      Matrix of exogenous variables
%
% by B. Kolb, Jan. 2015, based on code by J.P. LeSage (vare.m) 
%                        and A. Cesa Bianchi (VARmodel.m)


%% Housekeeping
if size(data,2) > size(data,1)
    data = data';
end

[T N] = size(data);

% Default determ
if ~exist('determ','var')
    determ = 1;
end

% Check if there are exogenous variables
if exist('varx','var')
    [Tx, Nx] = size(varx);
    % Check that data and varx are conformable
    if (Tx ~= T)
        error('Number of periods differs for data and exogenous variables!');
    end
else
    Nx = 0;
end

% Default nlagx
if ~exist('nlagx','var')
    nlagx = 0;
end

% Some useful statistics
ncoe  = N*nlags;
ncoex = Nx*(nlagx+1);
nobse = T - max(nlags,nlagx);
ncoet = ncoe + ncoex + determ;


%% Save some parameters and create data for VAR estimation
VARresults.nobs       = nobse;
VARresults.nvar       = N;
VARresults.nexog      = Nx;
VARresults.nlags      = nlags;
VARresults.nlags_exo  = nlagx;
VARresults.ncoeff     = ncoe;
VARresults.ncoeff_exo = ncoex;
VARresults.ncoeff_tot = ncoet;
VARresults.determ     = determ;

% Create independent vector and dependent regressor matrix of lags
[y, x] = VARform(data,nlags,determ);

% Create matrix of exogenous regressors (lags)
if ncoex > 0
    x_x  = VARform(varx,nlagx);
    if nlags == nlagx
        x = [x x_x];
    elseif nlags > nlagx
        x = [x x_x((nlags-nlagx)+1:end,:)];
    elseif nlags < nlagx
        y = y((nlagx - nlags)+1:end,:);
        x = [x((nlagx - nlags)+1:end,:) x_x];
    end
end


%% OLS estimation equation by equation
for j=1:N;
    yvec = y(:,j);
    OLS  = OLSreg(yvec,x);
    tout = tdis_prb(OLS.tstat,nobse-ncoe);
    % save OLS results
    eqnum = num2str(j);
    eval( ['VARresults.eq' eqnum '.beta  = OLS.beta;'] );
    eval( ['VARresults.eq' eqnum '.tstat = OLS.tstat;'] );
    eval( ['VARresults.eq' eqnum '.tprob = tout;'] );
    eval( ['VARresults.eq' eqnum '.resid = OLS.resid;'] );
    eval( ['VARresults.eq' eqnum '.yhat  = OLS.yhat;'] );
    eval( ['VARresults.eq' eqnum '.y     = yvec;'] );            
    eval( ['VARresults.eq' eqnum '.rsqr  = OLS.Rsqr;'] );
    eval( ['VARresults.eq' eqnum '.rbar  = OLS.Rbar;'] );
    eval( ['VARresults.eq' eqnum '.sige  = OLS.SSR;'] );
end


%% Compute the matrix of coefficients & VCV
Ft  = (x'*x)\(x'*y);
VCV = (y-x*Ft)'*(y-x*Ft)/(nobse-ncoet); % adjusted for # estim'd coeffs
VARresults.Ft    = Ft;
VARresults.sigma = VCV;
VARresults.resid = y - x*Ft;
VARresults.x     = x;
VARresults.y     = y;
if ncoex > 0
    VARresults.x_exo = x_x;
else
    VARresults.x_exo = 'NA';
end
