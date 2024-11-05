function [stateay, stateby, confidenceya, confidenceyb, bs_beta_dist, pval] = stateSVARLP_avg(hor, nlag, B, w, ...
    data,x,rpost,transformation,clevel,opt, ...
    y_tvar, I, rind, sind)

% NOTE: assumes no trend for the LP part of the irfs
% settings
h0 = 0; % hypothesis for the pvalues
bootstrap = 0; emp = 0; nstraps = 0; % to satisfy statelp func args
p = nlag;

[~, dsize] = size(data);
[T,K]= size(y_tvar);
rsize = length(rind);
if rsize == dsize
    r = dsize; % r = number of irfs response variables
else
    disp('error in number of IRF response variables')
end

% initialise matrices for later use

% POINT ESTIMATES
[lpa, lpb, ~, ~, ~, ~] = statelp(data,x,hor,rpost,transformation, clevel, opt, bootstrap, emp, p, nstraps);
[svara, svarb, ~, ~, ~, ~, beta] = stateSVAR(y_tvar,I,nlag,hor,rind,sind, B, 0,clevel);
stateay = w * lpa + (1-w) * svara;
stateby = w * lpb + (1-w) * svarb;

% BOOTSTRAPPING
% bootstrap settings
confidenceya = nan(2,r,hor);
confidenceyb = nan(2,r,hor);
bs_beta_dista = nan(B,r,hor);
bs_beta_distb = nan(B,r,hor);
bs_sea = nan(r,hor);
bs_seb = nan(r,hor);

% Generate residuals + make them dependent wild to bootstrap with
lagy_tvar = lagmatrix(y_tvar, [1:nlag]);
% (drop observations for all)
lagy_tvar = lagy_tvar(nlag+1:end,:);
ty_tvar = y_tvar(nlag+1:end,:);
tI = I(p+1:end,:);
TT = size(tI,1);
X = [ones(TT,1) lagy_tvar];

% (generate fitted and standard resid)
yhat_tvar = zeros(TT,K);
beta_a = squeeze(beta(:,:,1)); beta_b = squeeze(beta(:,:,2));
for i=1:TT
    yhat_tvar(i,:) = (X(i,:)*tI(i))*beta_a + (X(i,:)*(1-tI(i)))*beta_b;
end
res = ty_tvar - yhat_tvar;

% (convert to dependent wild)
eps = generate_dwb_resid(res, bw, dwb_settings);
% for i=1:K
%     eps(:,i) = generate_dwb_resid(res(:,i), bw, dwb_settings);
% end

eps = res;

for b=1:B
    y_bs = zeros(T,K); % we compute a T x K dependent variable via RFVAR DGP 
    y_bs(1:p,:) = y_tvar(1:p,:);
    res_samp = datasample(eps, TT, 'Replace', true);

    for j=p+1:TT
        ind = I(j);
        lagy_bs = [];
        for l=1:p
            lagy_bs = [lagy_bs, y_bs(j-l,:)];
        end
        x_bs = [1, lagy_bs];
        res = res_samp(j,:);
%         res = [eps(j,1) res_samp(j,2:3)];
        y_bs(j,:) = (x_bs*I(j))*beta_a + (x_bs*(1-I(j)))*beta_b + res;
    end
    % (state dependent LP irfs)
    data_temp = y_bs(:,rind);
    shock_temp = y_bs(:,sind);
    x_temp = lagmatrix(y_bs, [1:nlag]);
    
    shock_temp = shock_temp(nlag+1:end,:); % drop p observations
    data_temp = data_temp(nlag+1:end,:); % drop p observations
    x_temp = x_temp(nlag+1:end,:); % drop p observations (NaN)
    constant = ones(TT,1);
    x_temp = [constant, (1-tI).*constant, tI.*shock_temp, (1-tI).*shock_temp, repmat(tI,1,size(x_temp,2)).*x_temp, repmat((1-tI),1,size(x_temp,2)).*x_temp]; % assume no trend
    [lpa_temp, lpb_temp, ~, ~, ~, ~] = statelp(data_temp,x_temp,hor,rpost,transformation, clevel, opt, bootstrap, emp, nlag, nstraps);

    % (TVAR irfs)
    [svara_temp, svarb_temp, ~, ~, ~, ~] = stateSVAR(y_bs,I,nlag,hor,rind,sind, B, 0, clevel);
    
    % (store in irf distribution)
    bs_beta_dista(b,:,:) = w*lpa_temp + (1-w) * svara_temp;
    bs_beta_distb(b,:,:) = w*lpb_temp + (1-w) * svarb_temp;
end

for ri=1:r
    for hi=1:hor
        bs_sea(ri,hi) = std(bs_beta_dista(:,ri,hi));
        bs_seb(ri,hi) = std(bs_beta_distb(:,ri,hi));
    end
end

% Compute P values
pval = nan(2,r,hor);
taua = (stateay-h0)./bs_sea;
taub = (stateby-h0)./bs_seb;

for ri = 1:r
    for hi=1:hor
        % state A
        abs_taua = abs(taua(ri,hi)); % absolute value of analytical test statistic
        taua_bs = (bs_beta_dista(:,ri,hi) - (stateay(ri,hi)-h0))./bs_sea(ri,hi);
        abs_taua_bs = abs(taua_bs);
        pval(1,ri,hi) = (1/B) * sum(abs_taua_bs > abs_taua);

       % state B
       abs_taub = abs(taub(ri,hi));
       taub_bs = (bs_beta_distb(:,ri,hi)-(stateby(ri,hi)-h0))./bs_seb(ri,hi);
       abs_taub_bs = abs(taub_bs);
       pval(2,ri,hi) = (1/B) * sum(abs_taub_bs > abs_taub);
    end
end

% STORAGE
bs_beta_dist(1,:,:,:) = bs_beta_dista;
bs_beta_dist(2,:,:,:) = bs_beta_distb;

