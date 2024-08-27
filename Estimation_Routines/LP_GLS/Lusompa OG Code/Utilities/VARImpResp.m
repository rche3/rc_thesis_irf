
function [IRF, invA, IRF_opt] = VARImpResp(VAR,hor,ident,impact,S,shut)

% Compute impulse responses for a VAR model estimated with VARreg. Three
% identification schemes can be specified: zero short-run restrictions,
% zero long run restrictions, and sign restrictions
% 
% Usage:   [IRF IRF_opt] = VARImpResp(VAR,hor,ident,impact,S,shut)
%
% Inputs:  VAR:        Structure, as in output of VARreg function
%          hor:        Horizon over which impulse responses are generated
% (opt.)   ident:      Identification scheme; 'short_run' for zero short-
%                      run restrictions, 'long_run' for zero long-run 
%                      restrictions and 'sign' to generate one draw for 
%                      identification via sign restrictions 
%                      (default: 'short_run')
% (opt.)   impact:     Size of the shock; '0' for one standard deviation,
%                      '1' for unity (default: '0')
% (opt.)   S:          Random orthonormal matrix (if sign restriction
%                      identification is used) such that S*S'=I
% (opt.)   shut:       Index; shuts down any of the equations in the VAR
%                      for counterfactual analysis
%
% Outputs: IRF(h,v,s): Matrix of IRFs of variable v to shock s, of length h
%          Struct "IRF_opt" containing:
%  ~.hor:    Number of steps
%  ~.ident:  Identification chosen
%  ~.impact: Unit or 1 st dev
%  ~.invA:   Structural matrix such that invA*epsilon = u
%  ~.Fcomp:  Companion matrix
%  ~.maxEig: Max eigenvalue of the companion in absolute value
%
% by B. Kolb, Feb. 2015, based on VARir.m by A. Cesa Bianchi, see
%                        https://sites.google.com/site/ambropo/MatlabCodes

%% Check inputs
if ~exist('impact','var')
    impact = 0;
end

if ~exist('shut','var')
    shut = 0;
else
    if shut > VAR.nvar
        error('Input ''shut'' cannot be <= # variables.')
    end
end

if ~exist('ident','var')
    ident = 'short_run';
end

%% Retrieve parameters and preallocate IRFs
determ = VAR.determ;
nvar   = VAR.nvar;
nlags  = VAR.nlags;
Ft     = VAR.Ft;
sigma  = VAR.sigma;

IRF = zeros(hor,nvar,nvar);

%% Compute the companion matrix
F = Ft';
Fcomp = [F(:,1+determ:nvar*nlags+determ); eye(nvar*(nlags-1)) ...
    zeros(nvar*(nlags-1),nvar)];
Fcomp_eye = eye(size(Fcomp,1));
if shut~=0
    Fcomp(shut,:) = 0;
end

%% Compute the matrix invA containing the structural impulses
switch ident
    case 'short_run'
        [out, chol_flag] = chol(sigma);
        if chol_flag~=0; error('VCV is not positive definite'); end
        invA = out';
    case 'long_run'
        Finf_big = inv(eye(length(Fcomp))-Fcomp); % from the companion
        Finf     = Finf_big(1:nvar,1:nvar);
        D        = chol(Finf*sigma*Finf')'; % u2 has no LR effect on y1
        invA = Finf\D;
    case 'sign'
        [out, chol_flag] = chol(sigma);
        if chol_flag~=0; error('VCV is not positive definite'); end
        if ~exist('S','var'); error('No rotation matrix provided'); end
        invA = (S*out)';
end


%% Compute the impulse responses
for vv=1:nvar
    
    % Initialize impulse vector
    impulse  = zeros(nvar,1);
    % Set the size of the shock
    switch impact
        case 0 % one std. dev. shock
            impulse(vv,1) = 1;
        case 1 % unitary shock
            impulse(vv,1) = 1/invA(vv,vv);
        otherwise
            error('Input ''impact'' must be empty, 0 or 1.');
    end
    
    % Initialise impulse response matrix
    response = zeros(nvar, hor);
    % First period impulse response (to impulse vector)
    response(:,1) = invA*impulse;
    % Shut down the response of "shut"
    if shut~=0
        response(shut,1) = 0;
    end
    
    % Make impulse vector compatible with companion form of regressors
    impulse_big  = [(invA*impulse)' zeros(1, nvar*(nlags-1))]';
    
    % Recursive computation of impulse responses
    for hh = 2:hor
        Fcomp_eye = Fcomp_eye * Fcomp; % this is the multiplier Fcomp^n
        response_big   = Fcomp_eye * impulse_big;
        response(:,hh) = response_big(1:nvar);
    end
    IRF(:,:,vv) = response';
    Fcomp_eye = eye(size(Fcomp,1)); % reset for next iteration
end

%% Store used options
IRF_opt.hor    = hor;
IRF_opt.ident  = ident;
IRF_opt.impact = impact;
IRF_opt.invA   = invA;
IRF_opt.Fcomp  = Fcomp;
IRF_opt.maxEig = max(abs(eig(Fcomp)));
IRF_opt.shut   = shut;

