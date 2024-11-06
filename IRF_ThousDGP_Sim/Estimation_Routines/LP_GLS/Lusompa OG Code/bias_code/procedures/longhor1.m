function [vbias,betahat,betahat_adj,constant] = longhor1(yseries,xseries,first,last,nq,nk,narlag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: This is a routine.
% 
% Compute bias adjustment for long horizon regression nq periods ahead on lags 1 to nk of a regessor.  Invoked via
% 
%         execute longhor yseries xseries Wseries nWseries first last nq nk narlag vbias betahat betahat_adj
% 
%         execute longhor1 yseries xseries first last nq nk narlag vbias betahat betahat_adj
% 
% Passed by user in both routines:
% 
% yseries         univariate series for lhs variable
% xseries         univariate series for rhs variable(s).  Usually referred to as "x" in the description below.
% first           integer start date for regression
% last            integer end date for regression
% nq              integer horizon, nq>=0
% nk              integer number of lags of xseries to include on the rhs of the regression
% narlag          integer number of lags to include in estimating an AR model for xseries (needed
%                 to compute the bias)
% 
% Passed by user in longhor but not longhor1
% 
% Wseries         vector[series], in addition to xseries, used in computing moments related to xseries
% nWseries        the number of series in Wseries
% 
% Returned to user in both routines:
% 
% vbias           nkx1 vector, numerator of bias to order T (called "b" in the paper)
% bethat          nkx1 vector of regression coefficients
% bethat_adj      nkx1 vector, bias adjusted betahat, betahat_adj = betahat - (vbias/T)
% 
% The regression of interest is yseries, nq periods ahead, on lags 1 to nk of xseries:
% 
% (*)       yseries(t+nq) = const. + beta1*x(t-1) + ... + beta_nk*x(t-nk) + residual(t+nq)
% 
% where "x" is short for the parameter "xseries".  The regression is run with yseries dates running from first to last,
% i.e., with dates on x(t-1) running from first-nq-1 to last-nq-1.  The user specifies nq, first, last and nk.
% 
% Notes:
% 
% 1. The routine invokes procedures proc_vb_ma0, proc_vb_maq and proc_vbias, so these must be accesible from the program.
% 
% 2. The lhs variable may be a cumulated sum (yseries(t+nq) = dz(t)+dz(t+1)+...+dz(t+nq) (for a stationary series dz) or
% point in time.  If a cumulated sum, the user needs to appropriately construct yseries before invoking the routine.
% 
% 3. Obviously, (*) assumes that data on x(t) are available from observations first-nq-nk through last-nq-1.  But if
% data on x are also available from last-nq through last, the code will use that data as described in the next point.
% 
% 4. Within the routine, the moments necessary to compute bias rely in part from
% 
% *longhor1: a univariate AR in x(t) of order narlag.
% *longhor: a VAR in
% 
%         W(t)=(x(t),Wseries(1)(t), ..., Wseries(nWseries)(t))'
% 
% of order narlag.  W(t) is (nWseries+1) x 1.
% 
% The user specifies narlag; the routine estimates the autoregression or VAR. The code decides what sample to use in running the AR/VAR,
% via the following algorithm.  The sample runs from "firstx" and "lastx", defined as:
% 
% a. firstx: (1)If narlag<nk, the first date for x(t) on the lhs of the autoregression is firstx=first-nq-1, which is also the date of first
% observation on x(t-1) in regression (*).
% (2)If narlag>=nk, the date becomes firstx=first-nq+(narlag-nk).
% b. lastx: (1)If data on x are available in period "last", lastx=last. (Checked by seeing whether xseries(last)=%NA.)
% (2)If data on x(t) are not available in period "last", lastx=last-nq-1 -- that is, only use data on x also used in estimating
% regression (*).
% c. Obviously these dates can be adjusted, for example if one wants to experiment with different narlag but wants the sample to
% be the same for all narlag.
% 
% 5. The code assumes that the population disturbance is orthogonal to lagged x's (and Z's) and that the cumulants are zero.  The
% latter conditions rules out the disturbance being conditionally heteroskedastic or skewed.  See West (2016) on how to
% estimate bias in such a case.
% 
% 6. To better map into the notation of the paper:  The code defines
%                      yy(t) = yseries(t+nq)
% With this definition, the dating of the lhs variable yy(t) matches West (2016). Thus, for narlag>=1, nk>=1, and with W(t) defined
% above, regress
% 
%          yy(t) = const. + beta1*x(t-1) + ... + beta_nk*x(t-nk) + eta(t)
% 
%          x(t)  = const. + phi1*x(t-1) + ... + phi_narlag*x(t-narlag) + u(t) (longhor1)
% or
%          W(t) = const. + phi1*W(t-1) + ... + phi_narlag*W(t-narlag) + u(t) (longhor)
% 
% Note that, as in the paper, with this dating, eta(t) is actually realized in period t+nq.  To clarify: the parameters first and last
% passed by the user are the first and last observation dates on yseries (NOT yy).
% 
% 7. The code is largely undocumented.  But variable names map into the notation of West (2016) and other routines in a straightforward
% way.  The document "proc_vbias_basics" may be helpful in this respect.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the code proceeds in the following order
% (A)setup some parameters and variables
% (B)estimate AR of order narlag for x, construct companion form matrices phitwid and omegaUtwid
% (C)estimate the regression of yseries(t+nq) on x(t-1),...,x(t-nk)
% (D)compute EXeta and EetaZtwid, which are used in computing bias
% (E)comput bias, using external procedures proc_vb_ma0 (nq=0) or proc_vb_maq (nq>0)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (A)set up some parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% this construct firstx and lastx to prevent taking missing data
if narlag >= nk
    nZtwid = narlag;
firstx = first - nq + (narlag - nk);
else
    nZtwid = nk;
    firstx = first - nq - 1;
end

if (size(xseries,1)+(nq+1)) < last  %if there is missing data point, then
    lastx = last - nq - 1;
else
    lastx = last; 
end

phitwid = zeros(nZtwid,nZtwid);
omegaUtwid = zeros(nZtwid,nZtwid);
PX = eye(nk,nZtwid);
phinow = zeros(nZtwid,1);
EetaZtwid = zeros(nq+1,nZtwid);
EXeta = zeros(nk,1);
betahat = zeros(nk,1);
betahat_adj = zeros(nk,1);
vbias = zeros(nk,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (B)estimate AR of order nar for x, then use estimates to construct companion form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% to construct the relevant vector in VAR
con1 = ones(size(yseries,1),1); % constant term
nxlags = 1:narlag; % number of lags
nxonly1 = lagmatrix(xseries,nxlags); %lagged xseries
nX1 = [con1 nxonly1]; %regressors

% the AR regression
[beta_xseries,bint1,uresid] = regress(xseries(firstx:lastx,:),nX1(firstx:lastx,:));
clear bint1;

% phitwid, the companion form of phi, 
% is constructed by combining phinow (arranged coefficients) obtained from AR regression with identity matrix
for ii = 1:narlag
    phinow(ii,1) = beta_xseries(ii+1,1);
end

if nZtwid > 1
    phitwid(1,:) = phinow';
    phitwid(2:end,:)=eye(nZtwid-1,nZtwid);
else
    phitwid(1,1) = phinow(1,1);
end

% construct omegaUtwid
rss = uresid' * uresid;
omegaUtwid(1,1) = rss/size(uresid,1);
clear rss; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (C)estimate the regression of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% to construct the relevant vector in regression of interest
con2 = ones(size(yseries,1),1); % constant term
nxlags2 = nq+1:nq+nk;  % number of lags
nxonly2 = lagmatrix(xseries,nxlags2); % lagged xseries
nX2 = [con2 nxonly2]; % regressors

%regression of interest
[beta_YY,bint2,eta] = regress(yseries(first:last,:),nX2(first:last,:)); 

clear bint2;

Tnow = size(eta,1);

% compute betahat
for iii = 1:nk
    betahat(iii,1) = beta_YY(iii+1,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (D)compute EXeta and EetaZtwid, which include moments necessary to estimate bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% compute EetaZtwid
xsum=zeros(lastx+1,1); 
etaZ=zeros(nq+1,nZtwid);
%xsum(t) = x(t)+x(t+1)+...+x(t+nq-1)
for iq=0:nq
    if iq~=nq
        if lastx+iq<=lastx  % if taking value out of matrix bound, stop/add zero instead
            xsum(first-nq-nk+1:lastx)=xsum(first-nq-nk+1:lastx)+xseries((first-nq-nk+1)+iq:lastx+iq);
        else
            xsum(first-nq-nk+1:lastx)=xsum(first-nq-nk+1:lastx)+[xseries((first-nq-nk+1)+iq:lastx);zeros(iq,1)];
        end
    end
    for iz=1:nZtwid
        etaZ(iq+1,iz)=sum(eta.*xseries((first-nq)-iz+1+iq:(lastx-nq)-iz+1+iq));
        EetaZtwid(iq+1,iz)=etaZ(iq+1,iz)/Tnow;
    end
end

% compute etaxsum
etaxsum=zeros(lastx-(first-nq)+1,1);
for ik=1:nk
    etaxsum(ik)=sum(eta.*xsum(first-nq-ik+1:lastx-nq-ik+1)); %%%%
    EXeta(ik)=etaxsum(ik)/size(eta,1);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (E)compute bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute vbias
if nq == 0
    EetaZtwid0 = EetaZtwid(1,:);
    vbias = proc_vb_ma0(nZtwid,phitwid,omegaUtwid,nk,PX,EetaZtwid0);
else
    vbias = proc_vb_maq(nZtwid,phitwid,omegaUtwid,nk,PX,EetaZtwid,EXeta,nq);
end
% omegaUtwid and EetaZtwid affected by data size


% compute betahat_adj
betahat_adj = betahat - (vbias/Tnow);

% compute vbias/Tnow
vbias_over_Tnow=vbias/Tnow;

%This is Amaze adding a constant for the output
constant=beta_YY(1);

% Formatting
% disp('beta_hat_VAR in order of ');
% disp('constant_term beta_X_(t-1)...beta_X_(t-narlag)');
% disp([beta_xseries]);
%beta_xseries
%beta_YY
%vbias
%vbias_over_Tnow
%betahat
%betahat_adj
end